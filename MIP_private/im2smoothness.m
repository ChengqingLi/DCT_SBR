function s = im2smoothness(img, mode)
% s = im2smoothness(img, mode)
% This function returns the global smoothness (sum of differences between
% adjacent pixel values).
%
% Input arguments:
% img: an input image
% mode: 'h' = horizontal, 'v' = vertical
%       'hv' = both horizontal and vertical (default)
%
% Currently this function works with gray-scale images only and if the
% input is a color image it will be converted to gray-scale.
%
% Shujun Li @ University of Surrey, UK 2015
% Shujun Li: http://www.hooklee.com

s = -1; % Rerurn an invalid value by default.

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

if ischar(img)
    try
        img_filename = img;
        img = imread(img_filename);
    catch me
        fprintf('Failed to open the image file: %s!\n', me.message);
        return;
    end
end
if size(img,3)>1
    try
        img = rgb2gray(img);
        disp('The input is a color image. Converted to a gray-scale image!');
    catch me
        fprintf('Failed to convert the input color image into a gray-scale image: %s!\n', me.message);
        return;
    end
end

if ~exist('mode','var')
    mode = 'hv';
end

h = size(img,1);
w = size(img,2);

bHorizontal = true;
bVertical = true;
switch mode
    case 'h'
        bVertical = false;
    case 'w'
        bHorizontal = false;
end

if bVertical
    s = sum(sum(abs(diff(img))));
end

if bHorizontal
    s = s + sum(sum(abs(diff(img'))));
end
