function [outParams, time] = fmap_est_graphcut(imDataParams, algoParams)
% Claire Lin, Sept. 2020: add outputs to Hernando's graphcut for convergence comparison
% Description: Fat-water separation using regularized fieldmap formulation and graph cut solution. 
%
% Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.
% 
% Some properties:
%   - Image-space
%   - 2 species (water-fat)
%   - Complex-fitting
%   - Multi-peak fat (pre-calibrated)
%   - Single-R2*
%   - Independent water/fat phase
%   - Requires 3+ echoes at arbitrary echo times (some choices are much better than others! see NSA...)
%
% Input: structures imDataParams and algoParams
%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
%   - imDataParams.TE: echo times (in seconds)
%   - imDataParams.FieldStrength: (in Tesla)
%
%   - algoParams.species(ii).name = name of species ii (string)
%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%   Example
%      - algoParams.species(1).name = 'water' % Water
%      - algoParams.species(1).frequency = [0] 
%      - algoParams.species(1).relAmps = [1]   
%      - algoParams.species(2).name = 'fat' % Fat
%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
% 
%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%   - algoParams.range_r2star = [0 0]; % Range of R2* values
%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
%   - algoParams.range_fm = [-400 400]; % Range of field map values
%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
%   - algoParams.lambda = 0.05; % Regularization parameter
%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing) 
%
% Output: structure outParams
%   - outParams.species(ii).name: name of the species (taken from algoParams)
%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%   - outParams.r2starmap: R2* map (in s^{-1}, size [nx,ny])
%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%
%
% Author: Diego Hernando
% Date created: August 5, 2011
% Date last modified: November 10, 2011


DEBUG = 0;
fmiters = [];
time = [];


% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_graphcut( imDataParams,algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  outParams = [];
  return;
end

% Get data dimensions
[sx,sy,sz,C,N] = size(imDataParams.images);
% If more than one slice, pick central slice
if sz > 1
  disp('Multi-slice data: processing central slice');
  imDataParams.images = imDataParams.images(:,:,ceil(end/2),:,:);
end
% If more than one channel, coil combine
if C > 1
  disp('Multi-coil data: coil-combining');
  imDataParams.images = coilCombine(imDataParams.images);
end
  

% If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

% Check spatial subsampling option (speedup ~ quadratic SUBSAMPLE parameter)
SUBSAMPLE = algoParams.SUBSAMPLE;
if SUBSAMPLE > 1
  images0 = imDataParams.images;
  START = round(SUBSAMPLE/2);
  [sx,sy] = size(images0(:,:,1,1,1));
  allX = 1:sx;
  allY = 1:sy;
  subX = START:SUBSAMPLE:sx;
  subY = START:SUBSAMPLE:sy;
  imDataParams.images = images0(subX,subY,:,:,:);
end

% Regularization parameter
lambda=algoParams.lambda;

% Spatially-varying regularization.  The LMAP_POWER applies to the
% sqrt of the curvature of the residual, and LMAP_POWER=2 yields
% approximately uniform resolution.
LMAP_POWER = algoParams.LMAP_POWER;
  
  
% LMAP_EXTRA: Extra flexibility for including prior knowledge into
% regularization. For instance, it can be used to add more smoothing
% to noise regions (by adding, eg a constant LMAP_EXTRA), or even to
% add spatially-varying smoothing as a function of distance to
% isocenter...
LMAP_EXTRA = algoParams.LMAP_EXTRA;

% Finish off with some optimization transfer -- to remove discretization
DO_OT = algoParams.DO_OT;

% Let's get the residual. If it's not already in the params, compute it
try 
  % Grab the residual from the params structure
  residual = algoParams.residual;
catch
  % Check for uniform TE spacings
  dTE = diff(imDataParams.TE);
  
  try 
    TRY_PERIODIC_RESIDUAL = algoParams.TRY_PERIODIC_RESIDUAL;
  catch
    TRY_PERIODIC_RESIDUAL=1;
  end
  
  if TRY_PERIODIC_RESIDUAL==1 && sum(abs(dTE - dTE(1)))<1e-6 % If we have uniform TE spacing
    UNIFORM_TEs = 1;
  else
    UNIFORM_TEs = 0;
  end    


  if DEBUG == 1
    UNIFORM_TEs
  end

  % Compute the residual
  if UNIFORM_TEs == 1 % TEST DH* 090801
    % Find out the period in the residual (Assuming uniformly spaced samples)
    dt = imDataParams.TE(2)-imDataParams.TE(1);
    period = abs(1/dt);
    NUM_FMS_ORIG = algoParams.NUM_FMS;
    range = diff(algoParams.range_fm);
    params.NUM_FMS = ceil(algoParams.NUM_FMS/range*period);
    params.range_fm = [0 period*(1-1/(algoParams.NUM_FMS))];
    residual = computeResidual( imDataParams, algoParams );
    num_periods = ceil(range/period/2);
    algoParams.NUM_FMS = 2*num_periods*algoParams.NUM_FMS;
    residual = repmat(residual,[2*num_periods 1 1]);
    algoParams.range_fm = [-num_periods*period (num_periods*period-period/NUM_FMS_ORIG)];
  else
    % If not uniformly spaced TEs, get the residual for the whole range
      residual = computeResidual( imDataParams, algoParams );
  end

end

%save tempres.mat residual


% Setup the estimation, get the lambdamap,...
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
lmap = getQuadraticApprox( residual, dfm );  
lmap = (sqrt(lmap)).^LMAP_POWER;
lmap = lmap + mean(lmap(:))*LMAP_EXTRA;

% Initialize the field map indices
cur_ind = ceil(length(fms)/2)*ones(size(imDataParams.images(:,:,1,1,1)));

% This is the core of the algorithm
%fm = graphCutIterations(imDataParams,algoParams,residual,lmap,cur_ind );
[fm, fmiters1, time1] = graphcut(imDataParams,algoParams,residual,lmap,cur_ind );

% If we have subsampled (for speed), let's interpolate the field map
if SUBSAMPLE>1
  fmlowres = fm;
  [SUBX,SUBY] = meshgrid(subY(:),subX(:));
  [ALLX,ALLY] = meshgrid(allY(:),allX(:));
  fm = interp2(SUBX,SUBY,fmlowres,ALLX,ALLY,'*spline');
  lmap = interp2(SUBX,SUBY,lmap,ALLX,ALLY,'*spline');
  fm(isnan(fm)) = 0;
  imDataParams.images = images0;
end
algoParams.lmap = lmap;

% Now take the field map fm and get the rest of the estimates
for ka=1:size(imDataParams.images,6) 
  
  curParams = imDataParams;
  curParams.images = imDataParams.images(:,:,:,:,:,ka);

  
  if algoParams.range_r2star(2)>0
    % DH* 100422 use fine R2* discretization at this point 
    algoParams.NUM_R2STARS = round(algoParams.range_r2star(2)/2)+1; 
    r2starmap(:,:,ka) = estimateR2starGivenFieldmap( curParams, algoParams, fm );
  else
    r2starmap(:,:,ka) = zeros(size(fm));
  end
  
  if DO_OT ~= 1
    % If no Optimization Transfer, just get the water/fat images
    amps = decomposeGivenFieldMapAndDampings( curParams,algoParams, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
    waterimage = squeeze(amps(:,:,1,:));
    fatimage = squeeze(amps(:,:,2,:));
    w(:,:,:,ka) = waterimage;
    f(:,:,:,ka) = fatimage;
  end  
end

% If Optimization Transfer is requested, do it now
if DO_OT == 1
  imDataParams.fmGC = fm;

  %algoParams.OT_ITERS = 10;
  algoParams.fieldmap = fm;
  algoParams.r2starmap = r2starmap;
  algoParams.lambdamap = sqrt(lambda*lmap);
  
  [outParams,fmiters2, time2] = optimtransfer( imDataParams, algoParams  );      
  fm = outParams.fieldmap;
  
  fmiters = cat(3,fmiters1,fmiters2);
  time = cat(1,time1,time1(end)+time2);

  % Now re-estimate the R2* map and water/fat images
  if algoParams.range_r2star(2)>0
    algoParams.NUM_R2STARS = round(algoParams.range_r2star(2)/2)+1; 
    r2starmap(:,:,ka) = estimateR2starGivenFieldmap( curParams,algoParams, fm );
  else
    r2starmap(:,:,ka) = zeros(size(fm));
  end
  
  amps = decomposeGivenFieldMapAndDampings( curParams,algoParams, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
  waterimage = squeeze(amps(:,:,1,:));
  fatimage = squeeze(amps(:,:,2,:));
  w(:,:,:,ka) = waterimage;
  f(:,:,:,ka) = fatimage;
  
end
if isempty(fmiters)
    fmiters = fmiters1;
    time = time1;
end

% Put results in outParams structure
try
  outParams.species(1).name = algoParams.species(1).name;
  outParams.species(2).name = algoParams.species(2).name;
catch
  outParams.species(1).name = 'water';
  outParams.species(2).name = 'fat';
end  

outParams.species(1).amps = w;
outParams.species(2).amps = f;
outParams.r2starmap = r2starmap;
outParams.fieldmap = fm;
outParams.fms = fmiters;


end

%% 
function [fm, fmiters, time] = graphcut(imDataParams,algoParams,residual,lmap,cur_ind )

DEBUG = 0;
SMOOTH_NOSIGNAL = 1; % Whether to "homogenize" the lambdamap after
                     % some iterations, to get a smoother fieldmap in
                     % low-signal regions

STARTBIG = 1;

dkg = 15; % After dkg iterations, we may switch to a more homogeneous
          % regularization, to achieve more smoothness in noise-only
          % regions

DISPLAY_ITER = 0;


% Initialize some auxiliary variables
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:))*(imDataParams.FieldStrength)];
lambda = algoParams.lambda;
dt = imDataParams.TE(2)-imDataParams.TE(1);
period = 1/dt;
[sx,sy,N,C,num_acqs] = size(imDataParams.images);
fms = linspace(algoParams.range_fm(1),algoParams.range_fm(2),algoParams.NUM_FMS);
dfm = fms(2)-fms(1);
resoffset = [0:(sx*sy-1)]'*algoParams.NUM_FMS;
[masksignal,resLocalMinima,numMinimaPerVoxel] = findLocalMinima( residual, 0.06 );
numLocalMin = size(resLocalMinima,1);
stepoffset = [0:(sx*sy-1)]'*numLocalMin;
ercurrent = 1e10;
clear cur_ind2
% $$$ [A] = createExpansionGraphVARPRO_fast( residual, dfm, ones(sx,sy), algoParams.size_clique,ones(sx,sy), ones(sx,sy));

% cl: add time output
time = zeros(algoParams.NUM_ITERS+1,1);
fmiters = zeros(sx,sy,algoParams.NUM_ITERS+1);
tt = tic;

% Main graph-cut loop. At each iteration, a max-flow/min-cut problem is solved
if isfield(imDataParams,'finit') %cl: init w/ same winit as other methods, for conv compare
    fm = imDataParams.finit;
else
    fm = zeros(sx,sy);
end

fmiters(:,:,1) = fm;

for kg=1:algoParams.NUM_ITERS

  if kg == 1 & STARTBIG==1
    lambdamap = lambda*lmap;
    ercurrent = 1e10;
    prob_bigJump = 1;
  elseif (kg==dkg & SMOOTH_NOSIGNAL==1) | STARTBIG == 0
    lambdamap = lambda*lmap;
    ercurrent = 1e10;
    prob_bigJump = 0.5;
  end
  
  % Get the sign of the current jump
  cur_sign = (-1)^(kg);
  
  % Configure the current jump move
  if rand < prob_bigJump % If we're making a "jumpmin" move.
    cur_ind2(1,:,:) = cur_ind;
    repCurInd = repmat(cur_ind2,[numLocalMin,1,1]);
    if cur_sign>0
      stepLocator = (repCurInd+20/dfm>=resLocalMinima) & (resLocalMinima>0);
      stepLocator = squeeze(sum(stepLocator,1))+1 ;
      validStep = masksignal>0 & stepLocator<=numMinimaPerVoxel;
    else
      stepLocator = (repCurInd-20/dfm>resLocalMinima) & (resLocalMinima>0);
      stepLocator = squeeze(sum(stepLocator,1));
      validStep = masksignal>0 & stepLocator>=1;
    end
    nextValue = zeros(sx,sy);
    nextValue(validStep) = resLocalMinima(stepoffset(validStep) + stepLocator(validStep));
    cur_step = zeros(sx,sy);
    cur_step(validStep) = nextValue(validStep) - cur_ind(validStep); 
    
    if rand < 0.5
      nosignal_jump = cur_sign*round(abs(deltaF(2))/dfm);
    else
      nosignal_jump = cur_sign*abs(round((period - abs(deltaF(2)))/dfm));  
    end
    cur_step(~validStep) = nosignal_jump;
    
  else % If we're making a standard jump move.
    all_jump = cur_sign*ceil(abs(randn*3));
    cur_step = all_jump*ones(sx,sy);
    nextValue = cur_ind + cur_step;
    
    % DH* 100309: fix errors where fm goes beyond range
    if cur_sign>0
      cur_step(nextValue(:)>length(fms)) = length(fms) - cur_ind(nextValue(:)>length(fms));
    else
      cur_step(nextValue(:)<1) = 1 - cur_ind(nextValue(:)<1);
    end
    
    if DEBUG == 1
      all_jump
    end
  end
  
  
  % Create the huge adjacency matrix
  if norm(lambdamap,'fro') > 0
    [A] = createExpansionGraphVARPRO_fast( residual, dfm, lambdamap, algoParams.size_clique,cur_ind, cur_step);
  else 
    [A] = createExpansionGraphVARPRO_fast( residual, dfm, lambdamap, algoParams.size_clique,cur_ind, cur_step);
  end
  A(A<0)=0;
  
  % Solve the max-flow/min-cut problem (max_flow is a function of the matlabBGL library)
  [flowvalTS,cut_TS,RTS,FTS] = max_flow(A',size(A,1),1);
    
  % Take the output of max_flow and update the fieldmap estimate with
  % the best neighbor (among a set of 2^Q neighbors, where Q is the
  % number of voxels)
  cut1 = (cut_TS==-1); 
  cut1b = 0*cut1 + 1;
  cut1b(end) = 0;
  if DEBUG == 1
    compareGC = [sum(sum(A(cut1b==1, cut1b==0)))  ,  sum(sum(A(cut1==1, cut1==0)))]
  end
  
  if sum(sum(A(cut1b==1, cut1b==0))) <= sum(sum(A(cut1==1, cut1==0))) 
    cur_indST = cur_ind;

    if DEBUG == 1
      disp('Not taken');
    end
  else
    cut = reshape(cut1(2:end-1)==0,sx,sy);
    cur_indST = cur_ind + cur_step.*(cut);
    erST = sum(sum(residual(cur_indST(:)+resoffset(:)))) + dfm^2*lambdamap(1,1)*sum(sum(abs(diff(cur_indST,2,1)).^2)) + dfm^2*lambdamap(1,1)*sum(sum(abs(diff(cur_indST,2,2)).^2));
    if DEBUG == 1
      disp('Taken');
% $$$   pixelschanged = sum(sum(cut));
    end
  end
  
  % Update the field map (described as a map of indices in this function)
  prev_ind = cur_ind; 
  cur_ind = cur_indST;

% $$$   size(fms)
% $$$   mymax = max(cur_ind(:))
  
  cur_ind(cur_ind<1) = 1;
  cur_ind(cur_ind>length(fms)) = length(fms);

  fm = fms(cur_ind);
  
  if DISPLAY_ITER == 1
    % If we want to see the fieldmap at each step
    imagesc(fm,[-600 600]);axis off equal tight,colormap gray;colorbar;
    title(['Iteration: ' num2str(kg)],'FontSize',24);
    try 
      [kg fm(45,35) fms(cur_ind(45,35) + cur_step(45,35))]
    end
    drawnow;
  end
  fmiters(:,:,kg+1) = fm;
  time(kg+1) = toc(tt);
end

end


%===============================================================



function [outParams,fmiters,time] = optimtransfer( imDataParams, algoParams )

DEBUG_MODE = 0;

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_graphcut( imDataParams,algoParams );
if validParams==0
  disp(['Exiting -- data not processed']);
  outParams = [];
  return;
end

% If precession is clockwise (positive fat frequency) simply conjugate data
if imDataParams.PrecessionIsClockwise <= 0 
  imDataParams.images = conj(imDataParams.images);
  imDataParams.PrecessionIsClockwise = 1;
end

% Get recon parameters and images
gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency(1))*(imDataParams.FieldStrength)];
relAmps = algoParams.species(2).relAmps;
range_fm = algoParams.range_fm;
t = reshape(imDataParams.TE,[],1);
[sx,sy,sz,C,N,num_acqs] = size(imDataParams.images);
images = permute(imDataParams.images,[1 2 5 4 6 3]);

fm0 = algoParams.fieldmap;
r2starmap = algoParams.r2starmap;
NUM_ITERS = algoParams.OT_ITERS;
lambdamap = algoParams.lambdamap;


% Need to create my finite-difference matrix
if DEBUG_MODE==1
  tic
end
Dx = spalloc(sx*sy,sx*sy,2*sx*sy);
Dy = spalloc(sx*sy,sx*sy,2*sx*sy);
Dxy1 = spalloc(sx*sy,sx*sy,2*sx*sy);
Dxy2 = spalloc(sx*sy,sx*sy,2*sx*sy);
lmapvec = reshape(lambdamap,[],1);
for k=1:sx*sy
  if rem(k-1,sx) ~= 0
    Dx(k,k) = -min(lmapvec(k), lmapvec(k-1));
    Dx(k,k-1) = min(lmapvec(k), lmapvec(k-1));
  end
  if k-sx >= 1
    Dy(k,k) = -min(lmapvec(k), lmapvec(k-sx));
    Dy(k,k-sx) = min(lmapvec(k), lmapvec(k-sx));
  end
  if k-sx >= 1 & rem(k-1,sx) ~= 0
    Dxy1(k,k) = -min(lmapvec(k), lmapvec(k-sx-1));
    Dxy1(k,k-sx-1) = min(lmapvec(k), lmapvec(k-sx-1));
  end
  if k-sx >= 1 & rem(k+1,sx) ~= 0
    Dxy2(k,k) = -min(lmapvec(k), lmapvec(k-sx+1));
    Dxy2(k,k-sx+1) = min(lmapvec(k), lmapvec(k-sx+1));
  end
end
Dlambda = sqrt(2)*[Dx; Dy; Dxy1 ; Dxy2];
if DEBUG_MODE==1
  toc
end

% Need to get a upper bound on the Hessian of the data term 
% This should be easy, because data term is voxel-independent
Phi = getPhiMatrixMultipeak( deltaF,relAmps,t );
Gamma = Phi*inv(Phi'*Phi)*Phi';
[TM,TN] = meshgrid(reshape(t,[],1));
maxSignal = sum(sum(max(abs(images).^2,[],3),4),5);

if norm(r2starmap(:))==0
  B2 = 4*pi^2*sum(sum(abs(Gamma).*(TN-TM).^2));
  d2bound = B2*maxSignal;
else
  B2 = zeros(sx,sy);
  curGamma = zeros(sx,sy,N,N,num_acqs);
  for kx=1:sx
    for ky=1:sy
      for ka=1:num_acqs
        curPhi = Phi.*exp(-r2starmap(kx,ky,ka)*repmat(t(:),[1 2]));
        curGamma(kx,ky,:,:,ka) = curPhi*inv(curPhi'*curPhi)*curPhi';
        B2(kx,ky) = B2(kx,ky) + 4*pi^2*sum(sum(abs(squeeze(curGamma(kx,ky,:,:,ka))).*(TN-TM).^2));
      end
    end
  end
  d2bound = B2.*maxSignal;
end

D2data = spdiags(reshape(d2bound,[],1),0,sx*sy,sx*sy);
curHessTotal = Dlambda'*Dlambda + D2data;


fm = fm0;
% cl: add time output
time = zeros(NUM_ITERS,1);
fmiters = zeros(sx,sy,NUM_ITERS);
tt = tic;

for kit = 1:NUM_ITERS
  
  % Form the linear system... we need the gradient here
  d1 = zeros(sx,sy);
  for ka=1:num_acqs
    for kc=1:C
      for kn=1:N
        for km=1:N
          if norm(r2starmap(:))==0
            d1 = d1 - j*conj(images(:,:,kn,kc,ka)).*images(:,:,km,kc,ka)*Gamma(kn,km).*exp(j*2*pi*fm*(t(kn)-t(km)))*2*pi*(t(kn)-t(km));
          else
            d1 = d1 - j*conj(images(:,:,kn,kc,ka)).*images(:,:,km,kc,ka).*curGamma(:,:,kn,km,ka).*exp(j*2*pi*fm*(t(kn)-t(km)))*2*pi*(t(kn)-t(km));
          end
        end
      end
    end
  end
  d1 = real(d1); %Remove numerical errors that introduce imaginary component
  
  
  % Solve it
  curGradTotal = reshape(d1,[],1) + Dlambda'*Dlambda*fm(:);
  [fmstep,flag] = pcg(curHessTotal,-curGradTotal,1e-6,400);
  
  % Update field map
  fm = fm + reshape(fmstep,sx,sy);
  
  if DEBUG_MODE == 1
    curroughness = norm(Dlambda*fm(:))
    imagesc(reshape(fm,sx,sy),[-150 150]);drawnow
  end
  fmiters(:,:,kit) = fm;
  time(kit) = toc(tt);
end

% Now that we have the field map, estimate the water/fat images
for ka=1:num_acqs
  curParams = imDataParams;
  curParams.images = imDataParams.images(:,:,:,:,:,ka);
  amps = decomposeGivenFieldMapAndDampings( curParams,algoParams, fm,r2starmap(:,:,ka),r2starmap(:,:,ka) );
  w(:,:,:,ka) = squeeze(amps(:,:,1,:));
  f(:,:,:,ka) = squeeze(amps(:,:,2,:));
end


% Put results in outParams structure
try
  outParams.species(1).name = algoParams.species(1).name;
  outParams.species(2).name = algoParams.species(2).name;
catch
  outParams.species(1).name = 'water';
  outParams.species(2).name = 'fat';
end  
  
outParams.species(1).amps = w;
outParams.species(2).amps = f;
outParams.r2starmap = r2starmap;
outParams.fieldmap = fm;
end