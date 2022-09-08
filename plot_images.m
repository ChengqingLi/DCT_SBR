function plot_images(varargin)
% plot_images(img, img1, ...)
% plot_images(img, img1, ..., [rows cols])
% plot_images([rows cols], img, img1, ...)
% plot_images(img, title, img1, title1, ...)
% plot_images(img, title, img1, title1, ..., [rows cols])
% plot_images([rows cols], img, title, img1, title1, ...)
%
% Show one or more images and their PSNR/SSIM values.
% The first image is considered as the original image.
% PSNR and SSIM are calculated based on double, but histograms are shown by
% treaing all images as uint8 images to make the comparison easier.
% The images can have a dynamic range of either [0,1] or {0, ..., 255}.
%
% Shujun Li @ http://www.hooklee.com 2010-2015
% Junaid Jameel Ahmad @ University of Konstanz, Germany 2010

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

arg_index = 1;
img_index = 1;
imgs = cell(1,nargin);
names = cell(1,nargin);
Xmax = 255*ones(1,nargin);
rows = 0;
cols = 0;
while arg_index<=nargin
    if isnumeric(varargin{arg_index})
        if numel(varargin{arg_index})==2
            rows = varargin{arg_index}(1);
            cols = varargin{arg_index}(2);
            arg_index = arg_index + 1;
        else
            imgs{img_index} = varargin{arg_index};
            if max(imgs{img_index}(:))<=1
                Xmax(img_index) = 1;
            else
                Xmax(img_index) = 255;
            end
            arg_index = arg_index + 1;
            if (arg_index<=nargin && ischar(varargin{arg_index}))
                names{img_index} = varargin{arg_index};
                arg_index = arg_index + 1;
            elseif img_index==1
                names{img_index} = 'Original Image';
            else
                names{img_index} = ['Recovered Image #' num2str(img_index-1)];
            end
            img_index = img_index + 1;
        end
    elseif ischar(varargin{arg_index})
        names{img_index} = varargin{arg_index};
        arg_index = arg_index + 1;
        if arg_index<=nargin
            if isnumeric(varargin{arg_index})
                imgs{img_index} = varargin{arg_index};
                arg_index = arg_index + 1;
                img_index = img_index + 1;
            end
        end
    end
end
imgs(img_index:end) = [];
names(img_index:end) = [];
img_number = img_index - 1;

if img_number<1
    disp('At least one input argument should be an image!');
    return;
end

% If rows and cols are not defined, determine them automatically.
if (rows==0 || cols==0)
    if img_number<=3
        rows = 1;
        cols = img_number;
    elseif img_number<=6
        rows = 2;
        cols = ceil(img_number/2);
    elseif img_number<=9
        rows = 3;
        cols = ceil(img_number/3);
    else % if img_number>9
        rows = round(sqrt(img_number) *0.75);
        cols = ceil(img_number/rows);
    end
end

% Draw the first (original) image.
figure, subplot(rows, cols, 1);
imshow(imgs{1},[0 Xmax(1)]);
xlabel(names{1});
title(sprintf('PSNR = \\infty\nSSIM = 1\nsmoothness = %g', im2smoothness_border(imgs{1})));
% Draw the other image_number-1 images.
for i=2:img_number
    subplot(rows, cols, i);
    imshow(imgs{i},[0 Xmax(i)]);
    xlabel(names{i});
    title(sprintf('PSNR = %g\nSSIM = %g\nsmoothness = %g', psnr(imgs{i},imgs{1}), ssim(imgs{i},imgs{1}), im2smoothness_border(imgs{i})));
end

if exist('figureFullScreen.m', 'file')
    figureFullScreen(gcf);
end
