function ferror = compute_rmsd(argsError, x, varargin)
% Claire Lin, May 2020
% plot RMSD (in Hz) for paper
%
% input:
% x: image to compare with
% argsError [n,3] cell: n methods, each w/ name, time, x, each of length neval
%
% option:
% 'mask' - default true(size(x))
%
% output:
% ferror [n,1] cell: y coordinates of plots, each of length neval
%
% This routine computes RMSD over the "mask"
% (was entire image in the original paper).
%
% 2022-11-10 Jeff Fessler: add missing pi to Hz calculation, add mask

n = size(argsError,1);
ferror = cell(n,1);

arg.plot = 1;
arg.linestyle = {};
arg.step = 1; % 1 dot per s iters
arg.mask = true(size(x));
arg = vararg_pair(arg, varargin);

x = x(arg.mask);

if arg.plot == 1
    figure;
end
s = arg.step;
h = zeros(1,n);
for In = 1:n
    tc = argsError{In,2};
    nevalc = length(tc);
    xc = reshape(argsError{In,3},[],nevalc);
    xc = xc(arg.mask(:), :);
    ferror{In} = vecnorm(xc-x,2,1)/2/pi/sqrt(numel(x));
    if arg.plot == 1
        hold on
        h(In) = semilogy(tc(1),ferror{In}(1),'.-');
        if isempty(arg.linestyle)
            semilogy(tc,ferror{In},'-');
        else
            semilogy(tc,ferror{In},arg.linestyle{In})
        end
        semilogy(tc(1:s:end),ferror{In}(1:s:end),'.')
    end
end
legend(h, argsError(:,1))
end
