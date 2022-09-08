function [lb,ub] = dct2range(H, W, tmin, tmax)
% [lb,ub] = dct2range(H, W, tmin, tmax)
% This function returns the valid ranges of all coefficients of HxW DCT.
%
% Input arguments:
% H, W: height and width of the 2-D DCT (default: H = 8, W = H)
% tmin: the minimal pixel value (default: 0)
% tmax: the maximal pixel value (default: 255)
% For JPEG DCT coefficients please use tmin=-128 and tmax=127.
%
% Output arguments:
% lb: the lower bounds 
% ub: the upper bounds
%
% Note that lb(i)=-ub(i) always hold for i>0 due to the symmetry of the 2-D
% DCT.
%
% Shujun Li @ www.hooklee.com 2010-2015

if ~exist('H','var')
    H = 8;
end
if ~exist('W','var')
    W = H;
end
if ~exist('tmin','var')
    tmin = 0;
end
if ~exist('tmax','var')
    tmax = 255;
end

lb = zeros(H, W);
ub = zeros(H, W);
for u=1:H
    for v=1:W
        for k=1:H
            for l=1:W
                AA = dct2_A(u-1, v-1, k-1, l-1, H, W);
                if AA>0
                    lb(u,v) = lb(u,v) + AA*tmin;
                    ub(u,v) = ub(u,v) + AA*tmax;
                elseif AA<0
                    lb(u,v) = lb(u,v) + AA*tmax;
                    ub(u,v) = ub(u,v) + AA*tmin;
                end
            end
        end
    end
end
