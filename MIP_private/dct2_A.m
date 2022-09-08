function y = dct2_A(u, v, k, l, H, W)
% y = dct2_A(u, v, k, l, H, W)
% This function returns the A(u,v,k,l) term of 2-D HxW (heightXwidth) DCT.
% u, v, k, l should be 0-based indices.
% Default values: u = v = k = l = 0, H = 8, W = H
%
% Shujun Li @ www.hooklee.com 2010

y = [];

if ~exist('H','var')
    H = 8;
end
if ~exist('W','var')
    W = H;
end
if ~exist('u','var')
    u = 0;
end
if ~exist('v','var')
    v = 0;
end
if ~exist('k','var')
    k = 0;
end
if ~exist('l','var')
    l = 0;
end

y = dct_C(u,H) * dct_C(v,W) * cos(pi*u*(k+0.5)/H) * cos(pi*v*(l+0.5)/W);
