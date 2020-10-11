  function w = winit_pwls_pcg_ls(w, y, delta, smap, varargin)
% function w = winit_pwls_pcg_ls(w, y, delta, smap, varargin)
% Claire Lin, Sept. 2020
% PWLS on discretized ML estimate of winit
%
% min(w) weighted quadratic + R(w)
%| in
%|	w0  	[np 1]	rough estimate of winit
%|	y	    [np nc n]	n sets of measurements for nc coils
%|	delta	[1 n]	row vector of n echo time offsets
%|  smap    [np nc]   coil maps
%|
%| option
%|	stepper {'qs',# its}	monotonic line search parameters (def: {})
%|	niter			# of iterations (def: 30)
%|	maskR	[(np)]	logical reconstruction mask (required!)
%|	order			order of the finite diff matrix (def: 2)
%|	l2b             regularization parameter (2^) (def: -6)
%|	gammaType		CG direction: PR = Polak-Ribiere or FR = Fletcher-Reeves
%|	precon			Preconditioner: 'diag', 'chol', 'ichol'(def))
%|	reset			# of iterations before resetting direction (def: inf)
%|  df              delta f value in water-fat imaging (def: 0)
%|  relamp          relative amplitude in multipeak water-fat  (def: 1)
%|  tol             tolerance for ichol (def: 1e-3)
%|
%| out
%|	w   [np 1]	PWLS estimate of winit

arg.stepper = {};
arg.niter = 30;
arg.maskR = [];
arg.order = 2;
arg.l2b = -6;
arg.gammaType = 'PR';
arg.precon = 'ichol';
arg.reset = inf;
arg.df = 0;
arg.relamp = 1;
arg.tol = 1e-3; % for ichol
arg = vararg_pair(arg, varargin);

%% prepare variables and finite differencing matrices
w = double(w);
w0 = w;
y = double(y);

if isempty(arg.maskR)
	fail('maskR required')
end

% setup the monotonic line search parameters
if isempty(arg.stepper)
	arg.stepper = {'qs', 3}; % quad surr with this # of subiterations
end

% create the sparse regularization matrix
R = Reg1(arg.maskR, 'beta', 2^arg.l2b, 'order', arg.order, ...
	'distance_power', 2, 'type_diff', 'spmat', 'type_penal', ...
	'mat');
C = R.C;
if ~issparse(C)
	fail('CC = C^H * C is too slow if not sparse')
end
CC = C' * C;

% apply mask to data if necessary
if length(w) ~= sum(arg.maskR(:))
	w = w(arg.maskR,:);
	y = y(arg.maskR,:,:);
end

[np,nc,n] = size(y);

cost = zeros(arg.niter+1,1);

% check the data size matches the echo times
if n ~= size(delta,2), fail 'need delta to be [1 n]', end

%% calculate the magnitude and angles used in the data-fit curvatures
abss = abs(smap);
sjtotal = sum(abss.^2, 2); %[np,1]
angy = angle(y);
angs = angle(smap);
if arg.df
    Gamma = phiInv(arg.relamp, arg.df, delta); %[L,L]
end
set = 1; %coil = 1;
nset = cumsum(1:n-1);nset = nset(end);
wj_mag = zeros(np,nset,nc,nc);
d2 = zeros(1,nset);
ang2 = zeros(np,nset,nc,nc);
%%
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
%% compute |s_c s_d' y_dj' y_ci| /L/s * (tj - ti)^2
sjtotal(sjtotal==0) = 1; %avoid outside mask = 0
wj_mag = wj_mag./sjtotal;
if ~arg.df
    wj_mag = wj_mag/n;
end
dCC = diag(CC);
%% initialize projections and NCG variables
CCw = CC * w;
oldinprod = 0;
warned.dir = 0;
warned.step = 0;
%% begin iterations
fprintf('\n ********** ite_solve: NCG **********\n')
cost(1) = 0.5*sum(sum(wj_mag,[2:4]).*((w - w0).^2)) + norm(C*w,'fro');
fprintf(' ite: %d , cost: %f3\n', 0, cost(1)) 
%%
for iter=1:arg.niter
	%% compute the gradient of the cost function and curvatures
    hderiv = sum(wj_mag, [2:4]).* (w - w0);
    hcurv = sum(wj_mag, [2:4]);
    grad = hderiv + CCw;
	ngrad = -grad;
    % apply preconditioner
	switch arg.precon
	case 'diag'
		H = hcurv + dCC;
		npregrad = ngrad ./ H;
	case 'chol'
		%spparms('spumoni',2) % this will let you see if sparse Cholesky is used
		H = spdiag(hcurv) + CC;
        L = chol(H, 'lower');
		npregrad = L' \ (L \ ngrad);
	case 'ichol'
		H = spdiag(hcurv) + CC;
		alpha = max(max(sum(abs(H),2) ./ diag(H)) - 2,0); 
        L = ichol(H,struct('type','ict','droptol',1e-3,'diagcomp',alpha));
		npregrad = L' \ (L \ ngrad);
    otherwise
		npregrad = ngrad;
	end

	% compute CG direction
	newinprod = ngrad' * npregrad;
	newinprod = real(newinprod); %should be real, but just in case

	if oldinprod == 0 || mod(iter, arg.reset) == 0
		ddir = npregrad;
		gamma = 0;
		ngradO = ngrad;
	else
		if strcmp(arg.gammaType,'FR') % Fletcher-Reeves
			gamma = newinprod / oldinprod;
			ddir = npregrad + gamma * ddir;

		elseif strcmp(arg.gammaType,'PR') % Polack-Ribeir
			gamma = real((ngrad - ngradO)' * npregrad) / oldinprod;
			ngradO = ngrad;

			if (gamma < 0)
				printm('RESETTING GAMMA, iter=%d', iter)
				gamma = 0;
			end

			ddir = npregrad + gamma * ddir;

		end
	end
	oldinprod = newinprod;

	% check if correct descent direction
	if ddir' * grad > 0
		if ~warned.dir
			warned.dir = 1;
			warn 'wrong direction so resetting'
			printm('<ddir,grad>=%g, |ddir|=%g, |grad|=%g', ...
				ddir' * grad, norm(ddir), norm(grad))
		end
		% reset direction if not descending
		ddir = npregrad;
		oldinprod = 0;
	end

	% step size in search direction
	Cdir = C * ddir; % caution: can be a big array for 3D problems

	% compute the monotonic line search using quadratic cost
	CdCd = Cdir'*Cdir;
	CdCw = ddir'*CCw;
	step = 0;
    for is=1:arg.stepper{2}
        
        % compute the curvature and derivative for subsequent steps
        if step ~= 0
            hderiv = sum(wj_mag, [2:4]).* (w - w0);
            hcurv = sum(wj_mag, [2:4]);
        end

        % compute numer and denom of the Huber's algorithm based line search
        denom = (ddir.^2)' * hcurv + CdCd;
        numer = ddir' * hderiv + (CdCw + step * CdCd); %could be tiny bit faster?

        if denom == 0
            warn 'found exact solution??? step=0 now!?'
            step = 0;
        else
            % update line search
            step = step - numer / denom;
        end

    end
    % update the estimate and the finite differences of the estimate
	CCw = CCw + step * C' * Cdir;
	w = w + step * ddir;
    %
    cost(iter+1) = 0.5*sum(sum(wj_mag,[2:4]).*((w - w0).^2)) + norm(C*w,'fro');
    fprintf(' ite: %d , cost: %f3\n', iter, cost(iter+1)) 
end
end

function Gamma = phiInv(relamp,df,delta)
n = length(delta);
phi = [ones(n,1) sum(relamp.*exp(1i*delta(:)*df),2)]; %[n,2]
Gamma = phi*inv(phi'*phi)*phi';
end

