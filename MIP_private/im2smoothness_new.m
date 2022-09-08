function s = im2smoothness_new(img, intraFac, interFac, H, W)
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

if ~exist('H','var')
   H=8;
end

if ~exist('W','var')
   W=8;
end

h = size(img,1);
w = size(img,2);

img = double(img);

h_ind = 0; % For debugging purpose
s = 0;
% Separately handle horizontal and vertical pixel pairs (which help to
% save some computational time than using a single loop for both types
% of pixel pairs).
% Vertical pixel pairs.
x_ind2 = 1; % Initialize x-variable index2 (the larger index)
for j=1:w % Column by column
    for i=2:h % Row by row (without the first row)
        x_ind1 = x_ind2;
        x_ind2 = x_ind2 + 1;
        h_ind = h_ind + 1;
        if mod(i,H)==1
            s = s + interFac*abs(img(x_ind1)-img(x_ind2)); 
        else
            s = s + intraFac*abs(img(x_ind1)-img(x_ind2)); 
        end
    end
    x_ind2 = x_ind2 + 1; % For i=1 (the skipped first row)
end
% Horizontal pixel pairs.
x_ind2 = h; % Initialize x-variable index2 (the larger index)
for j=2:w % Column by column (without the first column)
    for i=1:h % Row by row
        x_ind2 = x_ind2 + 1;
        x_ind1 = x_ind2 - h;
        h_ind = h_ind + 1;
        if mod(j,W)==1
            s = s + interFac*abs(img(x_ind1)-img(x_ind2)); 
        else
            s = s + intraFac*abs(img(x_ind1)-img(x_ind2)); 
        end
    end
end
assert(h_ind == (2*w*h - w - h),'h_ind wrong in calculating smoothness!')