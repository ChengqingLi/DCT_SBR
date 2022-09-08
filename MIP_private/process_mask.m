function mask = process_mask(mask, h, w)
% mask = process_mask(mask, h, w)
% mask = process_mask(mask, [h w])
% This function process the input argument 'mask' (a 2-D logical array) to
% prepare for a number of functions.
% It 'mask' is a file name, it tries to load a mask from an image and then
% repeat it to of size hxw (image size) if needed.
% 
% Shujun Li @ University of Surrey, UK 2015
% Shujun Li: http://www.hooklee.com/

if ~exist('h','var')
    h = 8;
elseif numel(h)>=2
    if ~exist('w','var')
        w = h(2);
    end
    h = h(1);
end
if ~exist('w','var')
    w = h;
end

% If not given return an all-false matrix (all sign bits are unknown).
mask_default = false([h w]);

if ~exist('mask','var')
    mask = mask_default;
    return;
end
    
if ischar(mask)
    try
        % Load the mask image and convert to a logical matrix.
        mask = logical(imread(mask));
    catch
        mask = mask_default;
        disp('Failed to open the mask file!');
        return;
    end
end

H = size(mask,1);
W = size(mask,2);
if (mod(h,H)==0 && mod(w,W)==0)
    mask = repmat(mask, [h/H w/W]);
else
    disp('The size of the given mask does not match the image size! Return the loaded/given mask itself!');
end
