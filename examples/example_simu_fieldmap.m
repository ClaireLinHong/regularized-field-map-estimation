% Claire Lin, Sept. 2020
% Experiment A in "Efficient Regularized Field Map Estimation in 3D MRI",
% 2020 TCI
% 2022-11-07 Jeff Fessler small modifications to help debug Julia version
% ------------ Please setup MIRT before run -------------

%% add to path: data and functions
addpath('../data')
addpath('../functions')

%% variables in simulation
zp = 1:40; % choose subset of slices
nc = 4; % number of coils in simulation
% load data
load('input_object_40sl_3d_epi_snr40.mat')
mask = maskR(:,:,zp);
wtrue = double(in_obj.ztrue(:,:,zp)).*mask;
mag = double(in_obj.xtrue(:,:,zp)).*mask;

%% set parameters for data generation
p.etime = [0 2 10] * 1e-3; % echo times
p.true_thresh = 0.05; % threshold of the true object for determining reconstruction mask
p.yk_thresh = 0.1; % scale image
p.d_thresh = 0.1; % scale reg level
% initialize
ne = length(p.etime);
[nx,ny,nz] = size(mag);
yik = zeros(nx,ny,nz,nc,ne);

%% simulate sense map
smap = double(ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'nz', nz, 'ncoil', nc));
% normalize
tmp = sqrt(sum(abs((smap)).^2,4));
smap = div0(smap,tmp);

%% generate data
% multi-coil version
for kk = 1:ne
    yik(:,:,:,:,kk) = mag ...
        .* exp(1i * wtrue * (p.etime(kk) - p.etime(1))) ...
        .* smap;
end

% add complex Gaussian noise to image data
p.SNR = 24; % set noise level
p.SNR = inf; % todo
rng(0)
% compute the noise_std to get the desired SNR
image_power = 10*log10(sum(mag.^2,1:3)/(nx*ny*nz));
noise_power = image_power - p.SNR;
noise_std = sqrt(10^(noise_power/10));
noise_std = noise_std / 2; % because complex
% add the noise to the data
yikNN = yik;
yNse = noise_std * (randn(size(yik)) + 1i * randn(size(yik)));
yik = yik + yNse;
% compute the SNR
snr = zeros(1,ne);
for i = 1:ne
	tmp = norm(col(yikNN(:,:,:,:,i))) ./ norm(col(yNse(:,:,:,:,i)));
	snr(i) = mag2db(tmp);
end
pr snr
tmp = reshape(yik,[nx,ny,nz,nc,ne]);figure(1);im(tmp(:,:,:,1,1))

%% rescale yik
yik_sos = reshape(sum(yik.*reshape(conj(smap),[nx,ny,nz,nc]),4),[],ne); %coil combine
[yik_sos_scaled, scale] = ir_mri_field_map_reg_scale(yik_sos, p.etime, ...
		'fmax', p.yk_thresh, 'dmax', p.d_thresh);
yik_sos_scaled = reshape(yik_sos_scaled, [nx,ny,nz,ne]);


% 2d slice by slice
yik_scale = reshape(yik, [nx*ny, nz, nc, ne]) / scale;


%% Initialize fieldmap (winit):
% phase diff
printm(' -- Using simple initialization (no smoothing) and set of normalized first two images')
winit = angle( ...
    stackpick(yik_sos_scaled,2) .* conj(stackpick(yik_sos_scaled,1))) ...
    / (p.etime(2) - p.etime(1)); % difference of just first two sets
% set background pixels to mean of "good" pixels.
mag1 = abs(stackpick(yik_sos_scaled,1));
good = mag1 > p.yk_thresh * max(mag1(:));
winit(~good) = mean(winit(good));
% plot winit
winit_masked = winit(mask);
flim = [-70 70]; % brain range
figure(1); im(3, embed(winit_masked,mask)/2/pi, 'finit', flim), cbar('Hz')
clear mag1 good
drawnow


%% run QM-Huber / NCG for 3D, all slices
% initialize
yik_c = reshape(yik_scale, [nx*ny*nz,nc,ne]);
smap_c = reshape(smap, [nx*ny*nz,nc]);
l2b = -4;


%% 1. QM implementation
if ~isvar('wmap_qm')
 niter_qm = 250;
 [out,cost_qm,time_qm] = fmap_est_qm(winit(mask), yik_c(mask,:,:), p.etime,...
    smap_c(mask,:),'order', 1, 'l2b', l2b, 'dim',3, ...
    'niter', niter_qm,'maskR', mask,'hess','diag');
 wmap_qm = embed(out.ws,mask);
 figure(2);subplot(121);im(wmap_qm(:,:,:,end)/2/pi, 'qm', flim)
 subplot(122);semilogy(time_qm,cost_qm,'.-k')
 argsError_qm = {'QM', time_qm, wmap_qm};
end


%% 2. NCG implementation: no precon
if ~isvar('wmap_cg')
 niter_cg = 50;
%isave_cg = 0:niter_cg;
 [out,cost_cg,time_cg] = fmap_est_pcg_ls(winit(mask), yik_c(mask,:,:), ...
    p.etime, smap_c(mask,:), 'order', 1, 'l2b', l2b, ...
    'niter', niter_cg, 'maskR', mask, 'precon', 'none', 'gammaType', 'PR');
 wmap_cg = embed(out.ws, mask);
 fmap_cg = wmap_cg(:,:,:,end)/2/pi;
 figure(2); subplot(121); im('col', 8, fmap_cg, 'CG:I', flim); cbar
 subplot(122); semilogy(time_cg,cost_cg,'.-k')
 argsError_cg = {'NCG-MLS', time_cg, wmap_cg};
end


%% 3. NCG implementation: diag precon
if ~isvar('wmap_cg_d')
 niter_cg_d = 50;
 [out, cost_cg_d, time_cg_d] = fmap_est_pcg_ls(winit(mask), yik_c(mask,:,:), p.etime, ...
    smap_c(mask,:),'order', 1, 'l2b', l2b, ...
    'niter', niter_cg_d, 'maskR', mask, 'precon', 'diag', 'gammaType', 'PR');
 wmap_cg_d = embed(out.ws,mask);
 fmap_cg_d = wmap_cg_d(:,:,:,end)/2/pi;
 figure(2); subplot(121); im(fmap_cg_d, 'cg:diag', flim); cbar
 subplot(122); semilogy(time_cg_d, cost_cg_d, '.-k')
 argsError_cg_d = {'NCG-MLS-D', time_cg_d, wmap_cg_d};
end


%% 4. NCG implementation: ichol precon
if ~isvar('wmap_cg_i')
 niter_cg_i = 50;
 [out, cost_cg_i, time_cg_i] = fmap_est_pcg_ls(winit(mask), yik_c(mask,:,:),p.etime, ...
    smap_c(mask,:),'order', 1, 'l2b', l2b, ...
    'niter', niter_cg_i, 'maskR', mask, 'precon', 'ichol', 'gammaType', 'PR');
 wmap_cg_i = embed(out.ws,mask);
 fmap_cg_i = wmap_cg_i(:,:,:,end)/2/pi;
 figure(2); subplot(121); im(fmap_cg_i, 'cg:ic', flim)
 subplot(122); semilogy(time_cg_i, cost_cg_i, '.-k')
 argsError_cg_ichol = {'NCG-MLS-IC', time_cg_i, wmap_cg_i};
end


%% RMSE plots
argsError = {argsError_qm{:}; argsError_cg{:}; argsError_cg_d{:}; argsError_cg_ichol{:}};
error = compute_rmsd(argsError, wtrue, 'step', 10);
% labels
grid on
xlabel('Time (s)')
ylabel('RMSE (Hz)')
% color
ColorOdrDef = get(gca,'ColorOrder'); % 4x3 RGB array
ColorOdrCustom = [0 0 0;...
                0 1 0;...
                0 0 1;...
                1 0 0];
ColorOdrCustom = kron(ColorOdrCustom,[1;1;1]);
set(gca,'ColorOrder',ColorOdrCustom);
%% Image plots
zplot = 15:5:30;
scale = 1/2/pi;
clim = [-100,100];
wplot = wmap_cg_i(:,:,zplot,end);
figure;
subplot(511);
im(['notick'], 'row',1, [], [], mask(:,:,zplot) .* mag(:,:,zplot), [0,1], '|xtrue|'); cbar
subplot(512);
im(['notick'], 'row',1, [], [], winit(:,:,zplot)*scale, clim, 'f init'); cbar
subplot(513);
im(['notick'], 'row',1, [], [], wplot*scale, clim, 'f estimate'); cbar
subplot(514);
im(['notick'], 'row',1, [], [], mask(:,:,zplot).*wtrue(:,:,zplot)*scale, clim, 'true');cbar
subplot(515);
im(['notick'], 'row',1, [], [], mask(:,:,zplot).*(wtrue(:,:,zplot)-wplot)*scale, 5*[-10,10],'error');cbar
colormap gray
