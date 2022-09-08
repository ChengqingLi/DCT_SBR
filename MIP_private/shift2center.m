function y = shift2center(x, balanced_tails, Tmin, Tmax)
% y = shift2center(x, balanced_tails, Tmin, Tmax)
% 
% This function adjusts average amplitude of an input signal x towards the 
% center of the valid range [Tmin Tmax].
% y is the output signal after the amplitude adjustment.
% Note that the intensity adjustment is equaivalent to histogram shift, but
% we actually don't need to calculate the histogram at all.
% If balanced_tails=1, the histogram will be further shifted to balance the
% two tails (deadzones) around Tmin and Tmax.
% Default values: balanced_tails = 1, Tmin = 0, Tmax = 255.
%
% Copyright (C) 2010 Junaid Jameel Ahmad, Shujun Li (www.hooklee.com)

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

if ~exist('Tmax', 'var')
    Tmax = 255;
end
if ~exist('Tmin', 'var')
    Tmin = 0;
end
if Tmin>Tmax
    temp = Tmin;
    Tmin = Tmax;
    Tmax = temp;
end
if ~exist('balanced_tails', 'var')
    balanced_tails = 1;
end
Tmean = (Tmin + Tmax) / 2;

% Adjust the average intensity of the image to the center of the valid
% range [Tmin Tmax].
% This may cause overflow and/or underflow for some pixels.
x_mean = mean(x(:));
x = x + (Tmean - x_mean);

if balanced_tails % This can also be done when the histogram width is larger than Tmax-Tmin!
    % When the histogram width is smaller than the width of the valid range,
    % check if the intensity adjustment causes overflow or underflow.
    % If so, re-adjust the intensity to remove any overflow and underflow while
    % at the same time balance the tails of the histogram.
    y = x +  (Tmax-max(x(:))+Tmin-min(x(:)))/2;
else
    y = x;
end
