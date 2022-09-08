function y = idct2blockwise(x, H, W)
% y = idct2blockwise(x, H, W)
% y = idct2blockwise(x, [H W])
% This function performs blockwise 2-D WxH IDCT on the input array "x".
% Default values: H = 8, W = H
%
% Shujun Li @ www.hooklee.com 2010

y = [];

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

if ~exist('W','var')
    W = 8;
end
if ~exist('H','var')
    if numel(W)==1
        H = W;
    else
        H = W(2);
    end
end
H = H(1);

% Pad the image if needed to have an integral number of blocks.
if (mod(size(x,1),H)>0 || mod(size(x,2),W)>0)
    x = impad(x, W, H);
end

if H==W
    D = dctmtx(H);
end

x = double(x);
if exist('blockproc','file')==2
    if H==W
        fun = @(x) D'*x.data*D;
    else
        fun = @(x) idct2(x.data, [H W]);
    end
    y = blockproc(x, [H W], fun);
else
    if H==W
        fun = @(x) D'*x*D;
    else
        fun = @(x) idct2(x, [H W]);
    end
    y = blkproc(x, [H W], fun); %#ok<FPARK>
end
