 function [out, cost, times] = fmap_est_qm(w, y, delta, smap, varargin)
%function [out, cost, times] = fmap_est_qm(w, y, delta, smap, varargin)
%|
%| Phase unwrapping of multiple data sets (w denotes frequency)
%| using Huber's algorithm for quadratic surrogates
%| and sparse Cholesky factorization.
%| Can be 1d, 2d, etc., depending on R.
%| However, 3D problems may require too much memory for this to be practical!
%|
%| cost(w) = sum(i=0 to n-1) sum(j=0 to n-1)
%|		|yi*yj| (1 - cos(w*(dj-di) + \angle(yi) - \angle(yj)) + R(w)
%|
%| in
%|	w	[np 1]		initial estimate of field map w=omega
%|	y	[np nc n]	n sets of measurements for nc coils
%|	delta	[1 n]	row vector of n echo time offsets
%|                  note - the first delta is normally just 0.
%|  smap [np nc]    coil maps
%|
%| options
%|	niter			# of iterations (default 1)
%|	maskR	[(np)]	logical support mask (default full)
%|	order			order of the finite diff reg. (default 2)
%|	l2b             regularization parameter (2^)
%|  hess            'diag' or 'chol'
%|	dim             2 (2d) or 3 (3d) problem
%|  df              delta f value in water-fat imaging (def: 0)
%|  relamp          relative amplitude in multipeak water-fat  (def: 1)
%|
%| out
%|	out.ws	[np niter+1]	iterates
%|  out.xw / out.wf [np 1] water / fat images if arg.df~=0
%|	times	[niter+1 1]		time for each isave iteration (optional)
%|
%| 2020-01-31 Claire Lin, multicoil implementation

if nargin < 3, help(mfilename), error(mfilename), end

arg.niter = 1;
arg.maskR = [];
arg.order = 2;
arg.l2b = -6;
arg.hess = 'diag';
arg.dim = 3; % diag curvature for 2d / 3d
arg.df = 0;
arg.relamp = 1;
arg = vararg_pair(arg, varargin);

if isempty(arg.maskR)
	fail('Mask is required so that image dimensions are known!')
end
if ~islogical(arg.maskR)
	fail('Mask must be logical!')
end

w = double(w);
y = double(y);

% sparse finite-difference regularizing matrix (non-circulant)
R = Reg1(arg.maskR, 'beta', 2^arg.l2b, 'order', arg.order, ...
	'distance_power', 2, 'type_diff', 'spmat', 'type_penal', 'mat');
C = R.C;
if ~issparse(C)
	fail('CC = C^H * C is too slow if not sparse')
end
CC = C' * C;
if strcmp(arg.hess,'diag')
    if arg.order == 1
        dCC = 2^arg.l2b*4*arg.dim;
    elseif arg.order == 2
        dCC = 2^arg.l2b*16*arg.dim;
    else
        error('unknown order for C!')
    end
end

% apply mask to data if necessary
if numel(w) ~= sum(arg.maskR(:))
	w = w(arg.maskR);
	y = y(arg.maskR,:,:); % [np nc n]
end

[np,nc,n] = size(y);
times = zeros(arg.niter+1,1);
cost = zeros(arg.niter+1,1);

% check the data size matches the echo times
if n ~= size(delta, 2), fail 'need delta to be [1 n]', end

%% calculate the magnitude and angles used in the data-fit curvatures 
abss = abs(smap);
sjtotal = sum(abss.^2, 2); %[np,1]
angy = angle(y);
angs = angle(smap);
if arg.df
    Gamma = phiInv(arg.relamp,arg.df,delta); %[L,L]
end
set = 1; %coil = 1;
nset = cumsum(1:n-1);nset = nset(end);
wj_mag = zeros(np,nset,nc,nc);
d2 = zeros(1,nset);
ang2 = zeros(np,nset,nc,nc);
for j=1:n % for each pair of scans
    for i=1:n
        if i<j % only need one pair of differences
            d2(set) = delta(i) - delta(j);
            for c = 1:nc
                for d = 1:nc
                    wj_mag(:,set,c,d) = smap(:,c) .* conj(smap(:,d)) .*...
                        conj(y(:,c,i)) .* y(:,d,j); 
                    wj_mag(:,set,c,d) = abs(wj_mag(:,set,c,d));
                    % difference of the echo times and angles
                    ang2(:,set,c,d) = angs(:,c) - angs(:,d) + ...
                         angy(:,d,j) - angy(:,c,i);
                    if arg.df
                        wj_mag(:,set,c,d) = wj_mag(:,set,c,d)*abs(Gamma(i,j));
                        ang2(:,set,c,d) = ang2(:,set,c,d) + angle(Gamma(i,j));
                    end 
                end
            end
            set = set+1;
        end
    end
end
% compute |s_c s_d' y_dj' y_ci| /L/s * (tj - ti)^2
sjtotal(sjtotal==0) = 1; %cl: avoid outside mask 0 issue
wj_mag = wj_mag./sjtotal;
if ~arg.df
    wj_mag = wj_mag/n;
end
wm_deltaD = wj_mag .* d2;
wm_deltaD2 = wj_mag .* (d2.^2);
ang2(isnan(ang2))=0; %cl: avoid atan nan issue

% prepare output variables
out.ws = zeros(numel(w(:)), arg.niter+1);
out.ws(:,1) = w;

%% begin iterations
tt = tic; % start timing the iterations
fprintf('\n ********** ite_solve: QS Huber **********\n')
for iter = 1:arg.niter
    % compute num & denom contribution for each data-fit surrogate function
	% (numerator is the derivative of original cost function at current guess,
	% denominator is the curvatures from Funai)    
    [grad, denom, sm] = Adercurv(d2,ang2,wm_deltaD,wm_deltaD2,w);
    cost(iter) = sum(wj_mag.*(1-cos(sm)),'all') + norm(C*w,'fro');
    fprintf(' ite: %d , cost: %f3\n', iter-1, cost(iter)) 
    
	% add the regularization terms (and account for symmetric pairs with 2*)
	grad = 2*grad + (CC * w); % R.cgrad(R, w2);
    switch arg.hess
	case 'diag'
        H = 2*denom + dCC;
		tmp = grad ./ H;
    case 'chol'
        H = spdiag(2*denom) + CC;
        %tmp = H \ grad;
        L = chol(H, 'lower');
		tmp = L' \ (L \ grad);
    end
	w = w - tmp;

	% save iterates and time
	out.ws(:,iter+1) = w;
    times(iter+1) = toc(tt);

end

sm = w * d2 + ang2;
cost(iter+1) = sum(wj_mag.*(1-cos(sm)),'all') + norm(C*w,'fro');

fprintf(' ite: %d , cost: %f3\n', iter, cost(iter+1)) 
%output water & fat images
if arg.df
    x = decomp(w,arg.relamp,arg.df,delta,smap,y);
    out.xw = x(1,:).';
    out.xf = x(2,:).';
end
end 

function [hderiv, hcurv, sm] = Adercurv(d2,ang2,wm_deltaD,wm_deltaD2,w)
    % compute the data-fit derivatives and curvatures as in Funai paper
    sm = w * d2 + ang2;
    hderiv = sum(wm_deltaD .* sin(sm), [2:4]);

    srm = mod(sm + pi,2*pi) - pi;
    hcurv = sum(wm_deltaD2 .* ir_sinc_nopi(srm), [2:4]);
end

function Gamma = phiInv(relamp,df,delta)
n = length(delta);
phi = [ones(n,1) sum(relamp.*exp(1i*delta(:)*df),2)]; %[n,2]
Gamma = phi*inv(phi'*phi)*phi';
end

function x = decomp(w,relamp,df,delta,smap,y)
[np,nc,n] = size(y);
phi = [ones(n,1) sum(relamp.*exp(1i*delta(:)*df),2)]; %[n,2]
x = zeros(2,np);
for ip = 1:np
    B = phi.*col(exp(1i*w(ip)*delta(:))); %[n,2]
    B = kron(col(smap(ip,:)),B); %[n*nc,2]
    yc = permute(y(ip,:,:),[1,3,2]);
    x(:,ip) = B\yc(:);
end
end