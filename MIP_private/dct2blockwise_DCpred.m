function y = dct2blockwise_DCpred(z, DCpred_mode, H, W)
% y = dct2blockwise_DCpred(x, DCpred_mode, H, W)
% This function takes a matrix of DCT coefficients of an image and returns
% a new matrix where DC coefficients are differential encoded.
%
% Input arguments:
% DCpred_mode: mode of DC coefficient prediction (default: 3)
%       0: no prediction (return x directly)
%       1: predicted from previusly coded block in the same block-row
%       2: predicted from the previusly coded block in the same image
%          scanned using the raster order, as defined in JPEG standard)
%       3: predicted from one or two previusly coded blocks in the same
%          image scanned using the raster order, as defined in MPEG-1/2
%          standards)
% H, W: height and width of the 2-D DCT block (default: H = 8, W = H)
%
% Shujun Li @ University of Surrey, UK 2015
% Shujun Li: http://www.hooklee.com/

y = [];

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

h = size(z,1);
w = size(z,2);
if (mod(w,W)>0 || mod(h,H)>0)
    disp('The input DCT coefficient matrix does not match the block size!');
    return;
end

y = z;
switch(DCpred_mode)
    case 1
        for j=W+1:W:w
            for i=1:H:h
                y(i,j) = z(i,j) - z(i,j-W);
            end
        end
    case 2
        z_prev = 0;
        for i=1:H:h % Row by row
            for j=1:W:w % Column by column
                y(i,j) = z(i,j) - z_prev;
                z_prev = z(i,j);
            end
        end
    case 3
        % Process the first block-row.
        for j=W+1:W:w
            y(1,j) = z(1,j) - z(1,j-W);
        end
        % Process the first block-column.
        for i=H+1:H:h
            y(i,1) = z(i,1) - z(i-H,1);
        end
        % Process elements on other columns and rows.
        for j=W+1:W:w
            for i=H+1:H:h
                y(i,j) = z(i,j) - (z(i,j-W)+z(i-H,j))/2;
            end
        end
end

end
