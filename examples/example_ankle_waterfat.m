% Claire Lin, Sept. 2020
% Experiment D in "Efficient Regularized Field Map Estimation in 3D MRI",
% 2020 TCI
% ------------ Please setup MIRT before run -------------
%% add to path: data and functions 
addpath('./data')
addpath('./functions')
%% ankle data 
% ------------ Please download data from ISMRM fat-water data "USC/MultiChannel 3 Echo" -------------
% sensemap computed using ESPIRiT
load('Ankle_8ch.mat')
[nx,ny,nz,nc,ne] = size(data.images);
p.etime = data.TE;
df = 2*pi*440; %440 Hz for 3T
yik = data.images;
%% sensemap 
% ------------ Please estimate sensemap below (e.g. ESPIRiT) -------------
load('Ankle_sensemap.mat')
mag = sum(yik.*conj(sensemap),4);
mag0 = abs(stackpick(mag,1));
figure(1);im(mag0)
%% add complex Gaussian noise to image data
if 1
    p.SNR = 30; % set noise level 
    rng('default')
    rng(0)
    % compute the noise_std to get the desired SNR
    image_power = 10*log10(sum(mag0.^2,[1:3])/(nx*ny*nz));
    noise_power = image_power - p.SNR;
    noise_std = sqrt(10^(noise_power/10));
    noise_std = noise_std / 2; % because complex
    % add the noise to the data
    yikNN = yik;
    yNse = noise_std * (randn(size(yik)) + 1i * randn(size(yik)));
    yik = yik + yNse;
    % compute the SNR 
    SNR = zeros(ne,1);
    for i = 1:ne
        tmp = norm(col(yikNN(:,:,:,:,i))) ./ norm(col(yNse(:,:,:,:,i)));
        SNR(i) = mag2db(tmp);
    end
    pr SNR
    tmp = reshape(yik,[nx,ny,nz,nc,ne]);figure(1);im(tmp(:,:,:,1,1))
end
%% mask
p.true_thresh = 0.15;
maskErr = mag0 > p.true_thresh * max(mag0(:)); % uses the "true" magnitude
maskErr = logical(maskErr);
mask = zeros(nx,ny,nz);
for iz = 1:nz
    mask(:,:,iz) = bwconvhull(maskErr(:,:,iz), 'union');
end
mask = mask>0;
figure(1);im(mask)
%% Initialize fieldmap 1: discretization
printm(' -- Discretized ML estimation -- ')
yik_c = reshape(yik,[],nc,ne);
smap_c = reshape(sensemap,[],nc);
w0 = winit_water_fat(yik_c, p.etime, smap_c,'maskR',mask,'df',df,'nw',200);
w0 = embed(w0(:),mask);
wlim = 2*pi*[-200 100]; 
figure(2);subplot(131);im(w0,wlim)
% Initialize fieldmap 2: PWLS
printm(' -- PWLS on winit -- ')
y = reshape(yik,[],nc,ne);
w1 = winit_pwls_pcg_ls(w0(mask), yik_c(mask,:,:), p.etime,...
        smap_c(mask,:),'order', 1, 'l2b', -1,  ...
        'niter', 10,'maskR', mask,'precon','ichol',...
        'gammaType','PR','df',df);
w1 = embed(w1,mask);
figure(2);subplot(132);im(w1,wlim)
% Initialize fieldmap 3: set background pixels to mean of "good" pixels
winit = w1;
p.yk_thresh = 0.01;
good = (winit) > p.yk_thresh * max(yik(:));
winit(~good) = mean(winit(good));
winit = reshape(winit,[nx,ny,nz]);
figure(2); subplot(133);im(winit, wlim)
%% 1. QM implementation
niter_qm = 60;
[out,cost_qm,time_qm] = fmap_est_qm(winit(mask), yik_c(mask,:,:), p.etime,...
    smap_c(mask,:),'order', 1, 'l2b', -10, ...
    'niter', niter_qm,'maskR', mask,'hess','diag',...
    'dim',3,'df',df);
wmap_qm = embed(out.ws,mask);
%
figure(2);subplot(121);im(wmap_qm(:,:,:,end), 'qm', wlim)
subplot(122);semilogy(time_qm,cost_qm,'.-k')
argsError_qm = {'QM', time_qm, wmap_qm};
%% 2. NCG implementation: ichol precon
niter_cg = 10;
[out,cost_cg,time_cg] = fmap_est_pcg_ls(winit(mask), yik_c(mask,:,:),p.etime, ...
    smap_c(mask,:),'order', 1, 'l2b', -10, ...
    'niter', niter_cg,'maskR', mask,'precon','ichol',...
    'gammaType','PR','df',df);
wmap_cg = embed(out.ws,mask);
imwater_cg = embed(out.xw,mask);
imfat_cg = embed(out.xf,mask);
%
figure(2);subplot(121);im(wmap_cg(:,:,:,end), 'cg', wlim)
subplot(122);semilogy(time_cg,cost_cg,'.-k')
argsError_cg = {'NCG-MLS-IC', time_cg, wmap_cg};
%% RMSE plots
winf = wmap_cg(:,:,:,end);
argsError = {argsError_qm{:};argsError_cg{:}};
error = compute_rmsd(argsError, winf,'step',1);
% labels
grid on
axis([0,150,0,150])
xlabel('Time (s)')
ylabel('RMSD (Hz)')
% color
ColorOdrDef = get(gca,'ColorOrder'); % 4x3 RGB array
ColorOdrCustom = [0 0 0;...
                1 0 0];
ColorOdrCustom = kron(ColorOdrCustom,[1;1;1]);
set(gca,'ColorOrder',ColorOdrCustom);
%% Image plots
zplot = 1:4; 
figure;
% y and winit
subplot(611);im(['notick'], 'row',1, [], [], mask(:,:,zplot).*mag(:,:,zplot,1), [0,3],' ');cbar
clim = [-200,100];
subplot(612);im(['notick'], 'row',1, [], [], w0(:,:,zplot)/2/pi, clim,' ');cbar
subplot(613);im(['notick'], 'row',1, [], [], winit(:,:,zplot)/2/pi, clim,' ');cbar
% results
subplot(614);im(['notick'], 'row',1, [], [], wmap_cg(:,:,zplot,end)/2/pi, clim,' ');cbar
clim = [0,3];
subplot(615);im(['notick'], 'row',1, [], [], mask.*imwater_cg(:,:,zplot,1), clim,' ');cbar
subplot(616);im(['notick'], 'row',1, [], [], mask.*imfat_cg(:,:,zplot,1), clim,' ');cbar
colormap gray
