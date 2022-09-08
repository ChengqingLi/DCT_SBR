function s = im2smoothness_border_WA(img, H, W, wrapAround)
% s = im2smoothness_border_WA(img, H, W, wrapAround)
% This function returns the border smoothness (sum of differences between
% adjacent blocks).
%
% Input arguments:
% img: an input image
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
    H = 8;
end
if ~exist('W','var')
    W = H;
end
if ~exist('wrapAround','var')
    wrapAround = true;
end
h = size(img,1);
w = size(img,2);

img = double(img);

block_num_H = h/H;
block_num_W = w/W;
if wrapAround
    h_number = block_num_H*w + block_num_W*h; % number of h-variables = adjacent pixel pairs    
else
    h_number = (block_num_H-1)*w + (block_num_W-1)*h; % number of h-variables = adjacent pixel pairs
end
h_ind = 0; % For debugging purpose
s = 0;
% Separately handle horizontal and vertical pixel pairs (which help to
% save some computational time than using a single loop for both types
% of pixel pairs).
% Vertical pixel pairs.
x_ind2 = 1; % Initialize x-variable index2 (the larger index)
for j=1:w % Column by column
    for i=(H+1):H:h % Row by row (without the first row)
        x_ind2 = x_ind2 + H;
        x_ind1 = x_ind2 - 1;
        h_ind = h_ind + 1;
        s = s + abs(img(x_ind2)-img(x_ind1));
    end
    x_ind2 = x_ind2 + H; % For i=1 (the skipped first row)
end
if wrapAround
    s = s + 1/1.1*sum(abs(img(1,:)-img(h,:)));
    h_ind = h_ind + w;
end
Wh = W*h;
% Horizontal pixel pairs.
x_ind2 = Wh; % Initialize x-variable index2 (the larger index)
for j=(W+1):W:w % Column by column (without the first column)
    for i=1:h % Row by row
        x_ind2 = x_ind2 + 1;
        x_ind1 = x_ind2 - h;
        h_ind = h_ind + 1;
        s = s + abs(img(x_ind2)-img(x_ind1));
    end
    x_ind2 = x_ind2 + Wh - h;
end
if wrapAround
    s = s + 1/1.1*sum(abs(img(:,1)-img(:,w)));
    h_ind = h_ind + h;
end
assert(h_ind==h_number, 'h_ind wrong!');