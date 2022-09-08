function y = imscale(x, Tmax, Tmin)
% y = imscale(x, Tmax, Tmin)
% This function scale the input image x by mapping its maximum value and
% minimum value to Tmin and Tmax, respectively.
% The default values of Tmin and Tmax: Tmin = 0, Tmax = 255
% 
% Shujun Li @ www.hooklee.com 2011-2015

y = [];

if nargin<1
    disp('At least one input argument is needed!');
    return;
end
if ~exist('Tmax','var')
    Tmax = 255;
end
if ~exist('Tmin','var')
    Tmin = 0;
end

xmin = min(x(:));
xmax = max(x(:));

y = (x-xmin)/(xmax-xmin)*(Tmax-Tmin) + Tmin;

if (Tmax==255 && Tmin==0)
    y = uint8(y);
end
