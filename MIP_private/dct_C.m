function y = dct_C(u, N)
% y = dct_C(u, N)
% This function returns the C(u) constant in DCT.
% u should be 0-based index.
% Default values: u = 0, N = 8
%
% Shujun Li @ www.hooklee.com 2010

y = [];

if ~exist('N','var')
    N = 8;
end
if ~exist('u','var')
    u = 0;
end

if u==0
    y = sqrt(1/N);
else
    y = sqrt(2/N);
end
