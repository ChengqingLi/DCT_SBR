function z = dct2blockwise_DCrecov(y, DCpred_mode, H, W)
% y = dct2blockwise_DCrecov(x, DCpred_mode, H, W)
% This function takes a matrix of DCT coefficients of an image in which DC
% coefficients are differential encoded and returns the original DCT
% coefficients.
%
% Input arguments:
% DCpred_mode: mode of DC coefficient prediction (default: 3)
%              see dct2blockwise_DCpred for explanation of those methods
% H, W: height and width of the 2-D DCT block (default: H = 8, W = H)
%
% Shujun Li @ University of Surrey, UK 2015
% Shujun Li: http://www.hooklee.com/

z = [];

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

if ~exist('H','var')
    H = 8;
end
if ~exist('W','var')
    W = H;
end

if ~exist('DCpred_mode','var')
   DCpred_mode = 3;
end

h = size(y,1);
w = size(y,2);
if (mod(w,W)>0 || mod(h,H)>0)
    disp('The input DCT coefficient matrix does not match the block size!');
    return;
end

z = y;
switch(DCpred_mode)
    case 1
        for j=W+1:W:w
            for i=1:H:h
                z(i,j) = z(i,j) + z(i,j-W);
            end
        end
    case 2
        z_prev = 0;
        for i=1:H:h % Row by row
            for j=1:W:w % Column by column
                z(i,j) = z(i,j) + z_prev;
                z_prev = z(i,j);
            end
        end
    case 3
        % Process the first block-row.
        for j=W+1:W:w
            z(1,j) = z(1,j) + z(1,j-W);
        end
        % Process the first block-column.
        for i=H+1:H:h
            z(i,1) = z(i,1) + z(i-H,1);
        end
        % Process elements on other columns and rows.
        for j=W+1:W:w
            for i=H+1:H:h
                z(i,j) = z(i,j) + (z(i,j-W)+z(i-H,j))/2;
            end
        end
end

end
