 function w = winit_water_fat(y, delta, smap, varargin)
%function w = winit_water_fat(y, delta, smap, varargin) 
% Claire Lin, Sept. 2020
% Discretized ML estimate of fmap, voxel-wise
%| in
%|	y	[np nc n]	n sets of measurements for nc coils
%|	delta	[1 n]	row vector of n echo time offsets
%|  smap  [np nc]   sense maps
%| 
%| option
%| nw       number of discretization (def: 100)
%| relamp   relative amplitude in multipeak water-fat  (def: 1)
%| maskR	[(np)]	logical reconstruction mask (required!)
%| df               delta f value in water-fat imaging (def: 2*pi*440)

arg.nw = 100; 
arg.relamp = 1;
arg.df = 2*pi*440; % 220 Hz for 1.5T, 440 Hz for 3T
arg.maskR = [];
arg = vararg_pair(arg, varargin);
nw = arg.nw;
if isempty(arg.maskR)
	fail('maskR required')
end
y = double(y);
% apply mask to data if necessary
if size(y,1) ~= sum(arg.maskR(:))
	y = y(arg.maskR,:,:); % [np nc n]
    smap = smap(arg.maskR,:); % [np nc]
end
[np,nc,n] = size(y);
% check the data size matches the echo times
if n ~= size(delta,2), fail 'need delta to be [1 n]', end
%% calculate the magnitude and angles used in the data-fit curvatures
abss = abs(smap);
sjtotal = sum(abss.^2, 2); %[np,1]
angy = angle(y);
angs = angle(smap);
if arg.df
    Gamma = phiInv(arg.relamp,arg.df,delta); %[L,L]
end
set = 1; 
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
% compute |s_c s_d' y_dj' y_ci| /L/s 
sjtotal(sjtotal==0) = 1; 
wj_mag = wj_mag./sjtotal;
if ~arg.df
    wj_mag = wj_mag/n;
end
%%
% loop through increments
dfmax = max(abs(arg.df));
ws = linspace(-abs(dfmax/2),abs(dfmax/2),nw); 
cost = zeros(np,nw);
for iw = 1:nw
    if mod(iw,50) == 0
        fprintf(1,'Discretization # %d of %d\n', iw, nw);
    end
    w = ws(iw)*ones(np,1);
    sm = w * d2 + ang2;
    cost(:,iw) = sum(wj_mag.*(1-cos(sm)),[2:4]);
end
[~,idx] = min(cost,[],2);
w = ws(idx);
end

function Gamma = phiInv(relamp,df,delta)
n = length(delta);
phi = [ones(n,1) sum(relamp.*exp(1i*delta(:)*df),2)]; %[n,2]
Gamma = phi*inv(phi'*phi)*phi';
end