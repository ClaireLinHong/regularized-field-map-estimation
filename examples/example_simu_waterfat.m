% Claire Lin, Sept. 2020
% Experiment C in "Efficient Regularized Field Map Estimation in 3D MRI",
% 2020 TCI
% ------------ Please setup MIRT before run -------------
% ------------ Please download ISMRM fat-water toolbox to folder -------------
% (https://www.ismrm.org/workshops/FatWater12/data.htm)
%% add to path: data and functions 
addpath('./data')
addpath('./functions')
%% data 
% ------------ Please download data from ISMRM fat-water data "kellman_data" -------------
load('PKdata5.mat')
imDataParams = data; %[xyzcl]
imDataParams.images = double(data.images); 
%% algo parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0,0];%[0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 101; % Number of field map values to discretize
algoParams.SUBSAMPLE = 1;
algoParams.TRY_PERIODIC_RESIDUAL = 0;
%% add to path: golden-section search
path_lu = './fwtoolbox_v1_code/lu/';
addpath([path_lu]);
addpath([path_lu '/multiResSep']);
%% ground truth images
outParams2 = fw_3point_wm_goldSect(imDataParams, algoParams);
% plot images
imwater = outParams2.species(1).amps;
imfat = outParams2.species(2).amps;
ftrue = outParams2.fieldmap;
figure(1);
subplot(141);im(imwater)
subplot(142);im(imfat)
subplot(143);im(ftrue)
%% set parameters for data generation
p.te = data.TE; % echo times
p.SNR = 20; % set noise level 
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(data.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
[nx,ny,nz,nc,ne] = size(data.images);
% simu data
ytrue = zeros(nx,ny,nz,nc,ne);
for kk=1:ne
    ytrue(:,:,:,:,kk) = (imwater + imfat* ...
        sum(relAmps(:).*exp(1i*2*pi*deltaF(2:end).*p.te(kk)))) ...
        .*exp(1i*2*pi*ftrue*p.te(kk));
end
%% generate mask
[nx,ny,nz] = size(ftrue);
mag0 = abs(ytrue(:,:,:,:,1));
p.true_thresh = 0.05;
maskErr = mag0 > p.true_thresh * max(mag0(:)); % uses the "true" magnitude
maskErr = logical(maskErr);
mask = zeros(nx,ny,nz);
for iz = 1:nz
    mask(:,:,iz) = bwconvhull(maskErr(:,:,iz), 'union');
end
mask = mask>0;
%% simu data
yik = ytrue.*mask;
% add complex Gaussian noise to image data
if 1
    mag0 = abs(yik(:,:,:,:,1));
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
    SNR = zeros(ne,1);
    % compute the SNR 
    for i = 1:ne
        tmp = norm(col(yikNN(:,:,:,:,i))) ./ norm(col(yNse(:,:,:,:,i)));
        SNR(i) = mag2db(tmp);
    end
    pr SNR
end
% rescale
yik_long = reshape(yik,[],ne);
p.yk_thresh = 0.1; % scale image
p.d_thresh = 0.1; % scale reg level
[yik_scaled,scale] = ir_mri_field_map_reg_scale(yik_long, p.te, ...
		'fmax', p.yk_thresh, 'dmax', p.d_thresh);
yik_scaled = reshape(yik_scaled,[nx,ny,nz,nc,ne]);
imDataParams.images = yik_scaled.*mask; 
% 2d 
yik = reshape(yik_scaled,[nx*ny,nz,nc,ne]);
tmp = reshape(yik,[nx,ny,nz,nc,ne]);
%
figure(1);subplot(144);im(tmp(:,:,:,1,1))
%% Initialize fieldmap 1: discretization
smap = ones(nx,ny,nz,nc);
df = 2*pi*deltaF(2:end).'; %[1,ne]
wtrue = 2*pi*ftrue;
printm(' -- Discretized ML estimation -- ')
yik_c = reshape(yik,[],nc,ne);
smap_c = reshape(smap,[],nc);
w0 = winit_water_fat(yik_c, p.te, smap_c,'maskR',mask,...
    'relamp',relAmps,'df',df);
w0 = embed(w0(:),mask);
wlim = 2*pi*[-200 200]; 
figure(2);subplot(131);im(w0,wlim)
% Initialize fieldmap 2: PWLS
printm(' -- PWLS on winit -- ')
y = reshape(yik,[],nc,ne);
w1 = winit_pwls_pcg_ls(w0(mask), yik_c(mask,:,:), p.te,...
        smap_c(mask,:),'order', 2, 'l2b', -1,  ...
        'niter', 10,'maskR', mask,'precon','ichol',...
        'gammaType','PR','df',df,'relamp',relAmps);
w1 = embed(w1,mask);
figure(2);subplot(132);im(w1,wlim)
% Initialize fieldmap 3: set background pixels to mean of "good" pixels
winit = w1;
p.yk_thresh = 0.01;
good = winit > p.yk_thresh * max(yik(:));
winit(~good) = mean(winit(good));
% plot winit
winit = reshape(winit,[nx,ny,nz]);
figure(2); subplot(133);im(winit, wlim)
%% add to path: graphcut 
path_hernando = './fwtoolbox_v1_code/hernando/';
addpath([path_hernando 'common/']);
addpath([path_hernando 'descent/']);
addpath([path_hernando 'graphcut/']);
addpath([path_hernando 'mixed_fitting/']);
addpath([path_hernando 'create_synthetic/']);
addpath([path_hernando 'matlab-bgl-master/']);
%% 0. Graphcut implementation
algoParams.DO_OT = 0; % Optimal transfer
algoParams.OT_ITERS = 100; % number of optimal transfer iterations
algoParams.lambda = 2^(-7);% Regularization parameter
algoParams.LMAP_EXTRA = 0; % More smoothing for low-signal regions
algoParams.LMAP_POWER = 0; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.NUM_ITERS = 100; % Number of graph cut iterations
imDataParams.finit = winit/2/pi;
[outParams,time_gc] = fmap_est_graphcut(imDataParams, algoParams);
wmap_gc = 2*pi*outParams.fms.*mask;
argsError_gc = {'GC', time_gc, wmap_gc};
% plot image results
figure(3); 
imwater_gc = outParams.species(1).amps;
imfat_gc = outParams.species(2).amps;
subplot(231);im(ftrue)
subplot(232);im(imwater)
subplot(233);im(imfat)
subplot(234);im(wmap_gc(:,:,end))
subplot(235);im(imwater_gc)
subplot(236);im(imfat_gc)
% plot RMSE
argsError = {argsError_gc{:}};
error = compute_rmsd(argsError, wtrue);%,'linestyle',linestyle);
%% 1. QM implementation
niter_qm = 500;
[out,cost_qm,time_qm] = fmap_est_qm(winit(mask), yik_c(mask,:,:), p.te,...
    smap_c(mask,:),'order', 2, 'l2b', -7, ...
    'niter', niter_qm,'maskR', mask,'hess','diag',...
    'dim',3,'df',df,'relamp',relAmps);
wmap_qm = embed(out.ws,mask);
figure(2);subplot(121);im(wmap_qm(:,:,end), 'qm', wlim)
subplot(122);semilogy(time_qm,cost_qm,'.-k')
argsError_qm = {'QM', time_qm, wmap_qm};
%% 2. NCG implementation: ichol precon
niter_cg = 200;
[out,cost_cg,time_cg] = fmap_est_pcg_ls(winit(mask), yik_c(mask,:,:),p.te, ...
    smap_c(mask,:),'order', 2, 'l2b', -7, ...
    'niter', niter_cg,'maskR', mask,'precon','ichol',...
    'gammaType','PR','df',df,'relamp',relAmps);
wmap_cg = embed(out.ws,mask);
imwater_cg = embed(out.xw,mask);
imfat_cg = embed(out.xf,mask);
figure(2);subplot(121);im(wmap_cg(:,:,end), 'cg', wlim)
subplot(122);semilogy(time_cg,cost_cg,'.-k')
argsError_cg = {'NCG-MLS-IC', time_cg, wmap_cg};
%% RMSE plots
argsError = {argsError_qm{:};argsError_gc{:}; argsError_cg{:}};
error = compute_rmsd(argsError, wtrue,'step',20);
% labels
grid on
axis([0,100,0,60])
xlabel('Time (s)')
ylabel('RMSE (Hz)')
% color
ColorOdrDef = get(gca,'ColorOrder'); % 4x3 RGB array
ColorOdrCustom = [0 0 0;...
                0.4940 0.1840 0.5560;...
                1 0 0];
ColorOdrCustom = kron(ColorOdrCustom,[1;1;1]);
set(gca,'ColorOrder',ColorOdrCustom);
%% Field map plots
clim = [-100,200];
ie = 1;
figure;
subplot(131);im(['notick'], 'row',1, [], [], mask.*ytrue(:,:,:,:,ie), [0,150],' ');cbar
subplot(132);im(['notick'], 'row',1, [], [], w0/2/pi, clim,' ');cbar
subplot(133);im(['notick'], 'row',1, [], [], winit/2/pi, clim,' ');cbar
colormap gray
%% Image plots
figure;
% ground truth
subplot(531);im(['notick'], 'row',1, [], [], ftrue, clim,' ');cbar
subplot(532);im(['notick'], 'row',1, [], [], imwater,[0,100], ' ');cbar
subplot(533);im(['notick'], 'row',1, [], [], imfat,[0,100], ' ');cbar
% graph cut
subplot(534);im(['notick'], 'row',1, [], [], wmap_gc(:,:,end)/2/pi, clim,' ');cbar
subplot(535);im(['notick'], 'row',1, [], [], imwater_gc,[0,100], ' ');cbar
subplot(536);im(['notick'], 'row',1, [], [], imfat_gc,[0,100], ' ');cbar
subplot(537);im(['notick'], 'row',1, [], [], ftrue-wmap_gc(:,:,end)/2/pi, clim,' ');cbar
subplot(538);im(['notick'], 'row',1, [], [], imwater-imwater_gc,[0,100], ' ');cbar
subplot(539);im(['notick'], 'row',1, [], [], imfat-imfat_gc,[0,100], ' ');cbar
% ichol
subplot(5,3,10);im(['notick'], 'row',1, [], [], wmap_cg(:,:,end)/2/pi, clim,' ');cbar
subplot(5,3,11);im(['notick'], 'row',1, [], [], imwater_cg,[0,100], ' ');cbar
subplot(5,3,12);im(['notick'], 'row',1, [], [], imfat_cg,[0,100], ' ');cbar
subplot(5,3,13);im(['notick'], 'row',1, [], [], ftrue-wmap_cg(:,:,end)/2/pi, clim,' ');cbar
subplot(5,3,14);im(['notick'], 'row',1, [], [], imwater-imwater_cg,[0,100], ' ');cbar
subplot(5,3,15);im(['notick'], 'row',1, [], [], imfat-imfat_cg,[0,100], ' ');cbar
colormap gray