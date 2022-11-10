function error = compute_rmsd(argsError, x, varargin)
% Claire Lin, May 2020
% plot RMSD (in Hz) for paper
%
% input:
% x: image to compare with
% argsError [n,3] cell: n methods, each w/ name, time, x, each of length neval
% output:
% error [n,1] cell: y coordinates of plots, each of length neval
%
% This routine computes RMSD over the entire image, not just within the mask.
%
% 2022-11-20 Jeff Fessler: add missing pi to Hz calculation

n = size(argsError,1);
error = cell(n,1);
x = x(:);
arg.plot = 1;
arg.linestyle = {};
arg.step = 1; % 1 dot per s iters
arg = vararg_pair(arg, varargin);
if arg.plot == 1
    figure;
end
s = arg.step;
h = zeros(1,n);
for In = 1:n
    tc = argsError{In,2};
    nevalc = length(tc);
    xc = reshape(argsError{In,3},[],nevalc);
    error{In} = vecnorm(xc-x,2,1)/2/pi/sqrt(numel(x));
    if arg.plot == 1
        hold on;
        h(In) = semilogy(tc(1),error{In}(1),'.-'); 
        if isempty(arg.linestyle)
            semilogy(tc,error{In},'-');
        else
            semilogy(tc,error{In},arg.linestyle{In})
        end
        semilogy(tc(1:s:end),error{In}(1:s:end),'.')
    end
end
legend(h,argsError(:,1))
end
