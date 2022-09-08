function [imgs, Ts, imgs_dct2, PSNRs, SSIMs, BERs] = cplex_signbit(img_dct2, mask, DCpred_mode, DCTabs_min, relaxX, relaxZ, recovery_method, bDClevelshifted, H, W, Xmax)
% [imgs, Ts, imgs_dct2, PSNRs, SSIMs, BERs] = cplex_signbit(img_dct2, mask, DCpred_mode, DCTabs_min, relaxX, relaxZ, recovery_method, bDClevelshifted, H, W, Xmax)
% This function invokes IBM ILOG CPLEX's MATLAB Cplex class to solve the
% problem of recovering unknown sign bits of DCT coefficients of a
% DCT-transformed image "img_dct2" with a linear programming model.
%
% NOTE: Preliminary experiments showed C++ implementation is just slightly
% faster than the MATLAB implementation when the sign bits of the lowest 9
% DCT coefficient in each block are missing. More experiments are needed
% for more missing sign bits and for the deterministic times.
%
% Input arguments:
% img_dct2: the DCT-transformed image with unknown coefficients (can also
%           take an image or a file name as the input which will be
%           read/converted to its blockwise 8x8 2D DCT matrix)
% mask: a logical mask array showing which coefficients are unknown
%       (default: sign bits of all DC coefficients are unknown)
%       (can take a binary image or a filename as well)
% DCpred_mode: mode of DC coefficient prediction (default: 0)
%       0: no prediction
%       1, 2, 3: see dct2blockwise_DCpred
% DCTabs_min: threshold for removing DCT coefficients with unknown sign bits
%             but smaller amplitude (if amplitude <= DCTabs_min, the
%             coeffient will be set to 0 and not be recovered; default: 0)
%             NB: When DCTabs_min is too large it may create too strong
%             constraints so that a feasible solution does not exist!
% relaxX: true: x-variables are relaxed by xError for JPEG images. xError 
%             is calculated from the quantization errors of DCT
%             coefficients. false: x-variables are not relaxed. relaxX is
%             ignored if the input image is not JPEG. (default: true)
% relaxZ: true: z-variables are relaxed by quantization errors of JPEG
%             images. false: z-variables are not relaxed. relaxZ is ignored
%             if the input image is not JPEG. (default: true)
% recovery_method: the method for recovering the unknown sign bits
%                  0: all the methods (default)
%                  1: when the recovered DCT coefficient is 0, the result
%                     DCT coefficient is 0 as well
%                  2: when the recovered DCT coefficient is 0, the result
%                     DCT coefficient's sign is set to 1
%                  3: when the recovered DCT coefficient is 0, the result
%                     DCT coefficient's sign is set to 0
%                  4: when the recovered DCT coefficient is 0, the result
%                     DCT coefficient's sign is randomly assigned
%                  5: all signs are randomly guessed following random
%                     bounding, using a uniform distribution
%                  6: all signs are randomly guessed following random
%                     bounding, using a distribution learned from some test
%                     images
% bDClevelshifted: true = level shift DC coefficients by -round(xmax/2) as
%                  defined in JPEG standard (default)
%                  false = do not level shift DC coefficients (in this case
%                  DC coefficients' sign bits are always 1 so effectively
%                  known in all cases)
% H, W: height and width of the 2-D DCT block (default: H = 8, W = H)
% Xmax: maximal pixel value (default: 255)
%       The minimum pixel value Xmin is assumed to be always 0.
%
% Output arguments:
% imgs: a cell array holding a number of images recovered using different
%       methods
%       imgs.x0: the original image as input
%       imgs.x: the image recovered directly by CPLEX
%       imgs.x1: the image recovered from recovery method 1 (see the
%       explanation for the input argument recovery_method for more
%       details, the same for the next four fields)
%       imgs.x2: the image recovered from recovery method 2
%       imgs.x3: the image recovered from recovery method 3
%       imgs.x4: the image recovered from recovery method 4
%       imgs.x5: the image recovered from recovery method 5
%       imgs.x6: the image recovered from recovery method 6
%       imgs.x00: the image recovered from a naive method (all unknown
%                 signs are set to 0)
%       imgs.x01: the image recovered from a naive method (all unknown
%                 signs are set to 1)
%       imgs.x02: the image recovered from a naive method (all unknown
%                 signs are randomly assigned to 0 or 1)
% Ts: the time consumed by different part of the function
%     Ts.seconds_prep: time consumed by preprocessing input arguments in seconds
%     Ts.seconds_x0123: time consumed by creating initial condition(s)
%     Ts.seconds_model: time consumed by creating the CPLEX optimization model
%     Ts.seconds_cplex: time consumed by solving the CPLEX optimization model in seconds
%     Ts.ticks_cplex: deterministic time consumed by solving the CPLEX optimization model in ticks
%     Ts.seconds_solutions: time consumed by creating the final solutions in seconds
%     Ts.seconds_all: time consumed by the whole function in seconds
% imgs_dct2: DCT coefficients of some images (with differentially encoded
%            DC coefficients if DC prediction is involved)
%            imgs_dct2.y0: DCT coefficients of the input (original) image
%            imgs_dct2.y: DCT coefficients of img.x (the one recovered
%                         directly by CPLEX)
%            The above are useful for producing the figure on accuracy of
%            sign bit estimation using the optimization model.
% PSNRs, SSIMs: indicators of visual quality of recoverd images as measured
%               by PSNR and SSIM
% BERs: Bit Error Rates of different recovery methods
%       BER = Number of wrong bits / Total number of bits
% 
% If the input image is a JPEG image (when a filename is given as the first
% argument), then the function will handle the DCT coefficients according
% to the quantization matrix and the two's complement properly as the LP
% model will differ.
% 
% The mask can be defined to be of the same size as img_dct2, which allows
% different blocks to have different sets of unknown coefficients.
% If mask is of size HxW, it will be extended to be an array of the same
% size as the input image.
%
% Example: 
%     DCpred_mode=1, DCTabs_min=0, x-variables are not relaxed, z-variables are relaxed:
%     [imgs, Ts, imgs_dct2, PSNRs, SSIMs, BERs] = cplex_signbit('jpegs_quality\cameraman_100.jpg', 'masks\mask_9.pgm', 1, 0, false, true);
%
% Ruiyuan Lin, Shujun Li @ University of Surrey, UK 2010-2015
% Shujun Li: http://www.hooklee.com

t0 = tic; % Start time of the whole function.

imgs = [];
Ts = struct('seconds_prep',0,'seconds_model',0,'seconds_cplex',0,'ticks_cplex',0,'seconds_solutions',0,'seconds_x0123',0,'seconds_all',0);
imgs_dct2 = [];

img_filename = [];

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

xmin = 0;
if ~exist('xmax','var')
    Xmax = 255;
end
if (Xmax<=xmin)
    disp('xmax must be larger than xmin!');
    return;
end

if ~exist('bDClevelshifted','var')
    bDClevelshifted = true;
end
if bDClevelshifted
    DC_shift_level = round(Xmax/2); % e.g. 255 => 128 (JPEG standard)
    xmin = xmin - DC_shift_level;
else
    DC_shift_level = 0;
end
% Xmax is the maximum pixel value before level shift and xmax is the one
% after the level shift. Both are needed, so two variables.
% Xmin is always 0, so we don't need a separate variable for it.
xmax = Xmax - DC_shift_level;

if ~exist('H','var')
    H = 8;
end
if ~exist('W','var')
    W = H;
end
WH = W*H;

twosComp = false;

if ischar(img_dct2)
    try
        img_filename = img_dct2;
        img_dct2 = imread(img_filename);
        info = imfinfo(img_filename);
    catch me
        fprintf('Failed to open the image file: %s!\n', me.message);
        return;
    end
end

h = size(img_dct2,1);
w = size(img_dct2,2);
if (mod(w,W)>0 || mod(h,H)>0)
    disp('The input DCT-transformed image does not match the block size!');
    return;
end
wh = w*h;
wh2 = wh*2;
wh3 = wh*3;
wh4 = wh*4;

mask = process_mask(mask, h, w);
if ~bDClevelshifted
    if DCpred_mode == 0
        % If DC coefficients are not level shifted, then the DC coefficients'
        % sign bits are always 1 (known), so the mask should be adjusted.
        mask(1:H:h,1:W:w) = true;
    elseif DCpred_mode == 1
        mask(1:H:h,1) = true;
    elseif DCpred_mode == 2 || DCpred_mode == 3
        mask(1,1) = true;
    end
end
if all(mask(:)) % All-true matrix
    disp('No any sign bits are unknown, so nothing needs recovering!');
    return;
end

if ~exist('recovery_method','var')
   recovery_method = 0;
end
if ~exist('DCpred_mode','var')
   DCpred_mode = 0;
end
if ~exist('DCTabs_min','var')
   DCTabs_min = 0; 
end
if ~exist('relaxX','var')
    relaxX = true;
end
if ~exist('relaxZ','var')
    relaxZ = true;
end
if (isequal(class(img_dct2),'uint8') || max(img_dct2(:))<=1) % Input data is an image
    imgs.x0 = img_dct2;
    if exist('info','var') && strcmp(info.Format,'jpg')
        % If the input image is a JPEG image, handle the DCT coefficients
        % differently.
        twosComp = true;
        img_jpg = jpeg_read(img_filename);
        quan_z = img_jpg.coef_arrays{1,1};
        quan_y = dct2blockwise_DCpred(quan_z, DCpred_mode, H, W);
        quan_table = repmat(img_jpg. quant_tables{1,1}, [h/H w/W]);
        img_dct2 = quan_y .* quan_table;
    else
        % DC level shift should be applied to get the correct DC coefficients.
        img_dct2 = dct2blockwise(double(imgs.x0)-DC_shift_level, [H W]);
        % The DCT coefficients need adjusting if DC prediction is used.
        img_dct2 = dct2blockwise_DCpred(img_dct2, DCpred_mode, H, W);
    end
else
    % For img_dct2 directly input by users, we do not do DC level shift nor
    % DC prediction because we don't know if it is already processed.
    % Therefore, we assume it is DC level shifted and DC predicted.
    z = dct2blockwise_DCrecov(img_dct2, DCpred_mode, H, W);
    imgs.x0 = idctblockwise(z) + DC_shift_level;
end
imgs_dct2.y0 = img_dct2;
img_dct2_abs = abs(img_dct2);

Ts.seconds_prep = toc(t0);

t1 = tic;

% Pre-calculate dct2_A(u, v, k, l, H, W) for all possible values of
% (u,v,k,l) to avoid repeated computation later.
% (W*H)^2 cannot be too large, which is true for all existing image and 
% video coding standards (W = H = 8, 16 or 32).
data_filename = fullfile('data', ['dct2_As_' num2str(H) 'x' num2str(W) '.mat']);
% If the data is already in a file, simply load it which helps to save some
% time for re-calculating the static data.
if exist(data_filename,'file')
    load(data_filename);
end
if (exist('dct2_As','var') && isequal(size(dct2_As),[H W H W])) %#ok<NODEF>
    fprintf('Loaded needed static data of 2-D %dx%d DCT from ''%s''...\n', H, W, data_filename);
else
    fprintf('Pre-calculating and storing needed static data of 2-D %dx%d DCT...', H, W);
    dct2_As = zeros(H, W, H, W);
    for u=1:H
        for v=1:W
            for k=1:H
                for l=1:W
                    dct2_As(u, v, k, l) = dct2_A(u-1, v-1, k-1, l-1, H, W);
                end
            end
        end
    end
    % Save the data to avoid repeated calculation for future calls of the
    % same function.
    save(data_filename, 'dct2_As');
    disp('done.');
end

% ==================The 2-D DCT sign bits recovery model===================
% Variables:
% -- x(1)...x(wh) - pixel values
% -- y(1)...y(wh) - DCT coefficients (without DC prediction)
% -- QZ(1)...QZ(wh) - (quantized z).*(quantization table)
% -- z(1)...z(wh) - DCT coefficients (with DC prediction)
% -- h(x(i),x(j)) - absolute values between adjacent pixels, where x(i)
%                   and x(j) are adjacent pixel values
% Objective:
% -- minimize \sum_{x(i),x(j)} h(x(i),x(j))
% Constraints:
% -- y(i)=QZ(i) or y(i)=QZ(i)-QZ(j) or y(i)=QZ(i)-(QZ(j1)+QZ(j2))/2 for different
%    DC prediction modes
% -- QZ(i)-quan_table(i)/2<=z(i)<=QZ(i)+quan_table(i)/2
% -- x(k,l) = \sum_u\sum_v A(u,v,k,l)*z(u,v)
% -- x(i)-x(j) <= h(x(i),x(j))
% -- x(j)-x(i) <= h(x(i),x(j))
% -- y(u,v) = y*(u,v) for known y(u,v) whose value is y*(u,v)
% -- -|y(u,v)|<= y(u,v) <= |y(u,v)| if y(u,v)'s sign is unknown
% An equality x=a is represented by an inequality a<=x<=a.
% Lower and upper bounds:
% -- xmin <= x(i) <= xmax
% -- 0 <= h(x(i),x(j)) <= xmax
% y and z's lb and ub are determined implicitly via the DCT constraints.
% =========================================================================

% Error tolerance, a value smaller than 0.5 is sufficient to make sure
% rounding will work properly for x-variables (pixel values).
eps = 0.499;

cplex = Cplex(mfilename);
% cplex.Model.sense = 'minimize'; % minimization is the default.
% cplex.Param.lpmethod.Cur = 0; % 'auto' algorithm selection

fprintf('Setting the objective...');
h_number = wh2 - w - h; % number of h-variables = adjacent pixel pairs
var_number = wh4 + h_number;
% The variables: [{x(i)}_i {y(i)}_i {z(i)}_i {h(x(i),x(j)}_{x(i),x(j)}]
% Objective is 0*{x(i)}_i + 0*{y(i)}_i + 0*{z(i)}_i + {h(x(i),x(j)}_{x(i),x(j)}
cplex.Model.obj = [zeros(1,wh4) ones(1,h_number)];
disp('done.');

fprintf('Setting the lower and upper bounds of variables...');
% Set lb and ub for x-, z- and h-variables.
% Note that h-variables' lb-s are tight (0) so eps is not used.
% QZ-,z-variables will be uniquely determined by y-variables, so they don't
% need to have lb and ub set up.
if twosComp && relaxX
    % calculate the maximum error caused by JPEG DCT coefficient rounding
    xError = zeros(H,W);
    for i = 1:H
        for j = 1:W
            xError(i,j) = sum(sum(abs(dct2_As(:,:,i,j)) .* (img_jpg.quant_tables{1,1}/2))) + 1;
        end
    end
    xError = repmat(xError, [h/H w/W]);
    % calculate corresponding error tolerance for h-variables
    hError = zeros(1,h_number);
    h_ind = 0;
    x_ind2 = 1; % Initialize x-variable index2 (the larger index)
    for j=1:w % Column by column
        for i=2:h % Row by row (without the first row)
            x_ind1 = x_ind2;
            x_ind2 = x_ind2 + 1;
            h_ind = h_ind + 1;
            hError(h_ind)=xError(x_ind1)+xError(x_ind2);
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
            hError(h_ind)=xError(x_ind1)+xError(x_ind2);
        end
    end
    assert(h_ind==h_number,'h_ind wrong!');
    lb = [((xmin-eps)*ones(1,wh)-reshape(xError,1,wh)) -inf(1,wh3) zeros(1,h_number)];
    ub = [((xmax+eps)*ones(1,wh)+reshape(xError,1,wh)) inf(1,wh3) (Xmax+eps)*ones(1,h_number)+hError];
else
    lb = [(xmin-eps)*ones(1,wh) -inf(1,wh3) zeros(1,h_number)];
    ub = [(xmax+eps)*ones(1,wh) inf(1,wh3) (Xmax+eps)*ones(1,h_number)];
end
% Set lb and ub of y-variables (DCT coefficients with differentially
% encoded DC coefficients).
% Create two new variables to hold the lb and ub for later use (not in the
% LP model, but in the recovery methods).
lb_y = zeros(h,w);
ub_y = zeros(h,w);
for i=1:wh
    var_i = wh + i;
    if mask(i) % Known DCT coefficients
        lb_y(i) = img_dct2(i);
        ub_y(i) = img_dct2(i);
        lb(var_i) = img_dct2(i);
        ub(var_i) = img_dct2(i);
    else % DCT coefficients with unknown signs
        if twosComp && (mod(mod(i,h),H)==1 && mod(ceil(i/h),W)==1)% JPEG image, two's complement mode for DC coefficients.
            if quan_y(i)==0
                ub(var_i) = 0;
                lb(var_i) = 0;
              	lb_y(i) = 0;
                ub_y(i) = 0;
            else
                dc_dct_size = floor(log2(abs(quan_y(i)))) + 1;
                half_range = 2^(dc_dct_size - 1);
                full_range = half_range * 2;
                if quan_y(i)>0 % =0 is impossible
                    % Get the value of quan_y(i) without the sign bit.
                    exceptSignBit = quan_y(i) - half_range;
                    % Get the lb and ub of y-variables (for future use).
                    lb_y(i) = (exceptSignBit+1-full_range)*quan_table(i);
                    ub_y(i) = img_dct2(i);
                else % quan_y(i)<0
                    % Get the value of quan_y(i) without the sign bit.
                    exceptSignBit = quan_y(i) + full_range - 1;
                    % Get the lb and ub of y-variables (for future use).
                    lb_y(i) = img_dct2(i);
                    ub_y(i) = (exceptSignBit+half_range)*quan_table(i);
                end
                % Set the lb and ub of y-variables for the LP model.
                lb(var_i) = lb_y(i) - eps;
                ub(var_i) = ub_y(i) + eps;
            end
        else
            lb_y(i) = -img_dct2_abs(i);
            ub_y(i) = img_dct2_abs(i);
            lb(var_i) = lb_y(i) - eps;
            ub(var_i) = ub_y(i) + eps;
        end
        % Set the lb and ub of small y-variables to 0.
        % Note that the large value between the lower bound's absolute
        % value and the upper bound is used to compare with the threshold.
        % The two bounds can have different absolute values when two's
        % complement is used to encode differential DC coefficients.
        if max(-lb_y(i),ub_y(i))<DCTabs_min
            ub(var_i) = 0;
            lb(var_i) = 0;
        end
    end
end
cplex.Model.lb = lb;
cplex.Model.ub = ub;
imgs_dct2.lb_y = lb_y;
imgs_dct2.ub_y = ub_y;
disp('done.');

fprintf('Initializing the constraint matrix...');
con_number = wh3 + 2*h_number; % The number of constraints.
switch DCpred_mode
    % The number of non-zero elements in the constraints between y- and
    % z-variables depends on the DC pediction mode.
    case 0
        NZ_number0 = wh2;
    case {1,2}
        NZ_number0 = wh3;
    case 3
        NZ_number0 = wh4;
end
% The number of non-zero elements in the DCT constraints between x and z.
NZ_number1 = wh*(WH+1);
% The number of non-zero elements in the constraints between h and x.
NZ_number2 = 2*h_number*3;
% The number of non-zero elements in the constraints between h and x.
NZ_number3 = wh2;
% The maximal number of non-zero elements in the constraint matrix A.
NZ_number = NZ_number0 + NZ_number1 + NZ_number2 + NZ_number3;
% A vector of y-coordinates of all non-zero elements.
A_index_i = zeros(NZ_number,1);
% A vector of x-coordinates of all non-zero elements.
A_index_j = zeros(NZ_number,1);
% A vector of all non-zero elements.
A_value_ij = zeros(NZ_number,1);
if twosComp && relaxZ
    cplex.Model.lhs = [zeros(wh2,1); -quan_table(:)/2; -inf(2*h_number,1)];
    cplex.Model.rhs = [zeros(wh2,1); quan_table(:)/2; zeros(2*h_number,1)];
else
    cplex.Model.lhs = [zeros(wh3,1); -inf(2*h_number,1)];
    cplex.Model.rhs = zeros(con_number, 1);
end
% Initialize the number of constraints added.
NZ_index = 0;
disp('done.');

% Add the DCT constraint between x- and z-variables:
% x(k,l) = \sum_u\sum_v A(u,v,k,l)*z(u,v)
% This is actually represented here as:
% 0 <= -x(k,l) + \sum_u\sum_v A(u,v,k,l)*z(u,v) <= 0
% Note that here z(u,v) is the DCT coefficient after prediction.
fprintf('Setting the DCT constraints linking pixel values and DCT coefficients...');
for con_ind=1:wh
    i = mod(con_ind-1, h);
    j = floor((con_ind-1)/h);
    k = mod(i, H);
    l = mod(j, W);
    ii = i - k;
    jj = j - l;
    % Add the non-zero element corresponding to x(k,l).
    NZ_index = NZ_index + 1;
    A_index_i(NZ_index) = con_ind;
    A_index_j(NZ_index) = con_ind;
    A_value_ij(NZ_index) = -1;
    % Add the non-zero element corresponding to {y(u,v)}_{u,v}.
    A = dct2_As(:, :, k+1, l+1);
    for iii=1:H
        for jjj=1:W
            dct2_value = A(iii, jjj);
            if dct2_value~=0 % This helps save some computational time.
                NZ_index = NZ_index + 1;
                A_index_i(NZ_index) = con_ind;
                A_index_j(NZ_index) = wh3 + ((ii+iii)+(jj+jjj-1)*h);
                A_value_ij(NZ_index) = dct2_value;
            end
        end
    end
end
disp('done.');

% Add the constraints about DC prediction between y- and z-variables.
% y(i)=z(i) or y(i)=z(i)-z(j) or y(i)=z(i)-(z(j1)+z(j2))/2
% They are actually represented here as:
% 0 <= y(i) - z(i) <=0 (DCpred_mode = 0) or
% 0 <= y(i) - z(i) + z(j) <=0 (DCpred_mode = 1 and 2) or
% 0 <= y(i) - z(i) + 0.5*z(j1) + 0.5*z(j2) <=0 (DCpred_mode = 3)
fprintf('Setting the constraints about differentially encoded DC coefficients...');
Wh = W * h;
for j=1:w % Column by column
    for i=1:h % Row by row
        % Add the non-zero element corresponding to y(i).
        NZ_index = NZ_index + 1;
        con_ind = con_ind + 1;
        A_index_i(NZ_index) = con_ind;
        A_index_j(NZ_index) = con_ind; % = wh + sub2ind([h w],i,j)
        A_value_ij(NZ_index) = 1;
        % A variable used repeatedly below.
        wh_plus_con_ind = wh + con_ind; % = wh2 + sub2ind([h w],i,j)
        % Add the non-zero element corresponding to z(i).
        NZ_index = NZ_index + 1;
        A_index_i(NZ_index) = con_ind;
        A_index_j(NZ_index) = wh_plus_con_ind;
        A_value_ij(NZ_index) = -1;
        % Add the non-zero element corresponding to y(j) or y(j1) and
        % y(j2) according to the DC prediction mode.
        if (mod(i,H)==1 && mod(j,W)==1) % DC coefficients only, y(i)=z(i) for all AC
            switch DCpred_mode
                case 1
                    if j>1 % Skip the first of each block-row for which y(i)=z(i)
                        NZ_index = NZ_index + 1;
                        A_index_i(NZ_index) = con_ind;
                        A_index_j(NZ_index) = wh_plus_con_ind - Wh; % z(i,j-W)
                        A_value_ij(NZ_index) = 1;
                    end
                case 2
                    if j>1 % >=2nd DC coefficient of each block-row
                        NZ_index = NZ_index + 1;
                        A_index_i(NZ_index) = con_ind;
                        A_index_j(NZ_index) = wh_plus_con_ind - Wh; % z(i,j-W)
                        A_value_ij(NZ_index) = 1;
                    elseif i>1 % (and j==1): >=2nd DC coefficient of the first block-column (skip the 1st for which y(i)=z(i))
                        NZ_index = NZ_index + 1;
                        A_index_i(NZ_index) = con_ind;
                        A_index_j(NZ_index) = wh2 + (wh-Wh) + (i-H); % z(i-H,w-W+1)
                        A_value_ij(NZ_index) = 1;
                    end
                case 3
                    if (i==1 && j>1) % >=2nd DC coefficients of the first block-row
                        NZ_index = NZ_index + 1;
                        A_index_i(NZ_index) = con_ind;
                        A_index_j(NZ_index) = wh_plus_con_ind - Wh; % z(i,j-W)
                        A_value_ij(NZ_index) = 1;
                    elseif (j==1 && i>1) % >=2nd DC coefficients of the first block-column
                        NZ_index = NZ_index + 1;
                        A_index_i(NZ_index) = con_ind;
                        A_index_j(NZ_index) = wh_plus_con_ind - H; % z(i-H,j)
                        A_value_ij(NZ_index) = 1;
                    elseif (j>1 && i>1) % Other DC coefficients except the first one of the whole image for which y(i)=z(i))
                        NZ_index = NZ_index + 1;
                        A_index_i(NZ_index) = con_ind;
                        A_index_j(NZ_index) = wh_plus_con_ind - Wh; % z(i,j-W)
                        A_value_ij(NZ_index) = 0.5;
                        NZ_index = NZ_index + 1;
                        A_index_i(NZ_index) = con_ind;
                        A_index_j(NZ_index) = wh_plus_con_ind - H; % z(i-H,j)
                        A_value_ij(NZ_index) = 0.5;             
                    end
            end
        end
    end
end
disp('done.');

% Add the constraints between QZ- and z-variables.
% -quan_table(i)/2<=QZ(i)-z(i)<=quan_table(i)/2 for JPEG
% 0<=QZ(i)-z(i)<=0 for other image format
fprintf('Setting the constraints about quantization error...');
for j=1:w % Column by column
    for i=1:h % Row by row
        % Add the non-zero element corresponding to SB(i).
        NZ_index = NZ_index + 1;
        con_ind = con_ind + 1;
        A_index_i(NZ_index) = con_ind;
        A_index_j(NZ_index) = con_ind; % = wh2 + sub2ind([h w],i,j)
        A_value_ij(NZ_index) = 1;
        % A variable used repeatedly below.
        wh_plus_con_ind = wh + con_ind; % = wh3 + sub2ind([h w],i,j)
        % Add the non-zero element corresponding to z(i).
        NZ_index = NZ_index + 1;
        A_index_i(NZ_index) = con_ind;
        A_index_j(NZ_index) = wh_plus_con_ind;
        A_value_ij(NZ_index) = -1;
    end
end
disp('done.');

% Add constraints between x- and h-variables:
% x[ind1] - x[ind2] <= h[ind1,ind2]
% x[ind2] - x[ind1] <= h[ind1,ind2]
% which are equivalent to |x[ind1] - x[ind2]| <= h[ind1,ind2].
% They are actually represented here as:
% -inf < x[ind1] - x[ind2] - h[ind1,ind2] <= 0
% -inf < x[ind1] - x[ind2] - h[ind1,ind2] <= 0
% When \sum_{ind1,ind2}h[ind1,ind2] is minimized, we will have 
% h[ind1,ind2] = |x[ind1] - x[ind2]|.
fprintf('Setting the constraints between pixel values and the auxiliary h-variables...');
h_ind = wh4; % Initialize the index of the current h-variable
A_value_ij(NZ_index+(1:NZ_number2)) = repmat([1;-1;-1], [2*h_number 1]);
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
        var_indices = [x_ind1 x_ind2 h_ind; x_ind2 x_ind1 h_ind];
        for k=1:2
            con_ind = con_ind + 1;
            A_index_i(NZ_index+(1:3)) = con_ind;
            A_index_j(NZ_index+(1:3)) = var_indices(k,:);
            NZ_index = NZ_index + 3;
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
        var_indices = [x_ind1 x_ind2 h_ind; x_ind2 x_ind1 h_ind];
        for k=1:2
            con_ind = con_ind + 1;
            A_index_i(NZ_index+(1:3)) = con_ind;
            A_index_j(NZ_index+(1:3)) = var_indices(k,:);
            NZ_index = NZ_index + 3;
        end
    end
end
% Check if the number of variables, constraints, and non-zero elements
% in the constraint matrix are all valid.
assert(h_ind==var_number, 'The number of variables is wrong! Modelling error! Please check source code!');
assert(con_ind==con_number, 'The number of constraints is wrong! Modelling error! Please check source code!');
assert(NZ_index<=NZ_number, 'The number of none-zero elements in the constraint matrix is wrong! Modelling error! Please check source code!');
% Remove spare spaces reversed for more non-zero elements in the
% constraint matrix.
A_index_i(NZ_index+1:end) = [];
A_index_j(NZ_index+1:end) = [];
A_value_ij(NZ_index+1:end) = [];
% Create a sparse matrix to represent the constraint matrix.
cplex.Model.A = sparse(A_index_i, A_index_j, A_value_ij, con_number, var_number, NZ_index);
disp('done.');

Ts.seconds_model = toc(t1);

% Solve the optimization problem
disp('Solving the optimization problem...');
separator = repmat('=',1,80);
disp(separator);
cplex.solve();
disp(separator);

cplex_solution = cplex.Solution;

% Get real time in seconds for the core algorithm only
Ts.seconds_cplex = cplex_solution.time;
% Deterministic time in ticks (number of memory accesses, constant on
% different machines with the same computer architecture)
Ts.ticks_cplex = cplex_solution.dettime;

if ~isfield(cplex_solution,'x')
    disp('No solution is found!');
    return;
end

switch cplex_solution.method
    case 0
        cplex_method = 'None';
    case 1
        cplex_method = 'Primal simplex';
    case 2
        cplex_method = 'Dual simplex';
    case 4
        cplex_method = 'Barrier optimizer (no crossover)';
    case 11
        cplex_method = 'Feasopt';
    case 12
        cplex_method = 'Mixed integer optimizer';
    otherwise
        cplex_method = 'Others';
end
% The following is information for debugging.
imgs.x = cplex_solution.x(1:wh);
imgs.x = reshape(imgs.x, [h w]);
objvals = [cplex_solution.objval sum(cplex_solution.x(wh4+1:end)) im2smoothness(imgs.x)];
fprintf('Objective value of the returned solution: %g (Method: %s)\n', objvals(1), cplex_method);
fprintf('Objective value calculated from the retuned h-variables: %g\n', objvals(2));
fprintf('Objective value calculated from the retuned x-variables: %g\n\n', objvals(3));
assert(abs(objvals(1)-objvals(2))<eps && abs(objvals(2)-objvals(3))<eps, 'Objective values do not macth! CPLEX must have been running in an unstable mode! Results are not reliable!');

t1 = tic;

% Get the mask for DCT coefficients with unknown signs.
mask_neg = ~mask;
mask_neg_num = sum(mask_neg(:));
mask_neg_GtDCTabsMin = mask_neg & (max(-lb_y,ub_y) >= DCTabs_min);
mask_neg_GtDCTabsMin_num = sum(sum(mask_neg_GtDCTabsMin));

disp('Calculating the final solution(s)...');
% Get the directly recovered image by CPLEX.
imgs.x = imgs.x + DC_shift_level;
imgs.x = im_postprocess(imgs.x, Xmax);
% Recovered DCT coefficients (with differentially encoded DC coefficients)
imgs_dct2.y = cplex_solution.x(wh+1:wh2);
imgs_dct2.y = reshape(imgs_dct2.y, [h w]);
BERs.x = sum( imgs_dct2.y(mask_neg_GtDCTabsMin)~=img_dct2(mask_neg_GtDCTabsMin) )/mask_neg_GtDCTabsMin_num;

% Mask of DC coefficients with unknown signs and 0-value
mask_neg_Z = (mask_neg & imgs_dct2.y==0);
mask_neg_Z_num = sum(mask_neg_Z(:));
fprintf('The number of unknown sign bits that cannot be directly estimated: %d (%g%%).\n', mask_neg_Z_num, mask_neg_Z_num/mask_neg_num*100);
mask_neg_NZ = (mask_neg & imgs_dct2.y~=0);
% The variable y will be used repeatedly for all recovery methods.
y = img_dct2;
y(mask_neg_NZ) = (1 + sign(imgs_dct2.y(mask_neg_NZ)))/2 .* ub_y(mask_neg_NZ) + ...
                 (1 - sign(imgs_dct2.y(mask_neg_NZ)))/2 .*lb_y(mask_neg_NZ);
% Method 1: For DCT coefficients whose sign bits are still unknown
% (recovered DCT coefficients are 0), set those DCT coefficients to 0.
if (recovery_method==0 || recovery_method==1)
    if (mask_neg_Z_num>0)
        y(mask_neg_Z) = 0;
    end
    BERs.x1 = sum( y(mask_neg_GtDCTabsMin)~=img_dct2(mask_neg_GtDCTabsMin) )/mask_neg_GtDCTabsMin_num;
    imgs.x1 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
end
% Method 2: For DCT coefficients whose sign bits are still unknown
% (recovered DCT coefficients are 0), set their sign bits to 1
% (non-negative).
% If mask_neg_Z is empty, then imgs.x2 will be identical with imgs.x1, so
% no need to do this. The same for imgs.x3 and imgs.x4.
if (recovery_method==0 || recovery_method==2)
    if (mask_neg_Z_num>0)
        y(mask_neg_Z) = ub_y(mask_neg_Z);
        BERs.x2 = sum( y(mask_neg_GtDCTabsMin)~=img_dct2(mask_neg_GtDCTabsMin) )/mask_neg_GtDCTabsMin_num;
        imgs.x2 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
    else
        disp('All unknown sign bits can be directly estimated, so skip Recovery Method 2!');
    end
end
% Method 3: For DCT coefficients whose sign bits are still unknown
% (recovered DCT coefficients are 0), set their sign bits to 0
% (non-positive).
if (recovery_method==0 || recovery_method==3)
    if (mask_neg_Z_num>0)
        y(mask_neg_Z) = lb_y(mask_neg_Z);
        BERs.x3 = sum( y(mask_neg_GtDCTabsMin)~=img_dct2(mask_neg_GtDCTabsMin) )/mask_neg_GtDCTabsMin_num;
        imgs.x3 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
    else
        disp('All unknown sign bits can be directly estimated, so skip Recovery Method 3!');
    end
end
% Method 4: For DCT coefficients whose sign bits are still unknown
% (recovered DCT coefficients are 0), randomly guess their sign bits
% (50-50).
if (recovery_method==0 || recovery_method==4)
    if (mask_neg_Z_num>0)
        randNum = randi([0 1],mask_neg_Z_num,1);
        y(mask_neg_Z) = randNum .* ub_y(mask_neg_Z) + (1 - randNum) .* lb_y(mask_neg_Z);
        BERs.x4 = sum( y(mask_neg_GtDCTabsMin)~=img_dct2(mask_neg_GtDCTabsMin) )/mask_neg_GtDCTabsMin_num;
        imgs.x4 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
    else
        disp('All unknown sign bits can be directly estimated, so skip Recovery Method 4!');
    end
end
% Method 5: For all DCT coefficients, use random bounding with an
% assumed uniform distribution to decide their sign bits.
if (recovery_method==0 || recovery_method==5)
    % Calculate p1=Prob(SB=1) of each DCT coefficient to be the the
    % relative distance between the recovered value and the upper bound.
    % Alternatively, p0=Prob(SB=-1) is the relative distance between
    % the recovered value and the lower bound.
    % Assume uniform distribution of imgs_dct2.y on [lb_y ub_y].
    p1 = (imgs_dct2.y - lb_y)./ (ub_y - lb_y);
    randNum = (rand(mask_neg_num,1)<=p1(mask_neg));
    y(mask_neg) = randNum .* ub_y(mask_neg) + (1-randNum) .* lb_y(mask_neg);
    BERs.x5 = sum(y(mask_neg)~=img_dct2(mask_neg))/mask_neg_num;
    imgs.x5 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
end
% Method 6: For all DCT coefficients, use random bounding with a
% distribution learned from some test images to decide their sign bits.
% This method requires the distribution to be given, which will be
% done indirectly via a function d0_to_p1() mapping a relative distance
% between the recovered value and the lower bound (d0) to the p1 value.
% If the function does not exist, x6 will not be produced.
% The d0_to_p1() function should be able to handle a matrix element wise.
if ((recovery_method==0 || recovery_method==6) && exist('d0_to_p1.m','file'))
    d0 = (imgs_dct2.y - lb_y)./ (ub_y - lb_y);
    % d0_to_p1() function may behave differently for different parameters
    % including DCpred_mode, DC_shift_level, DCTabs_min and twosComp.
    p1 = d0_to_p1(d0, DCpred_mode, DC_shift_level, DCTabs_min);
    randNum = (rand(mask_neg_num,1)<=p1(mask_neg));
    y(mask_neg) = randNum .* ub_y(mask_neg) + (1-randNum) .* lb_y(mask_neg);
    BERs.x6 = sum(y(mask_neg)~=img_dct2(mask_neg))/mask_neg_num;
    imgs.x6 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
end
disp('done.');

Ts.seconds_solutions = toc(t1);


t1 = tic;
% Recover the image using a naive method (assigning all unknown sign bits to 1).
y = img_dct2;
y(mask_neg) = ub_y(mask_neg);
BERs.x01 = sum(y(mask_neg)~=img_dct2(mask_neg))/mask_neg_num;
imgs.x01 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
% Recover the image using a naive method (assigning all unknown sign bits to 0).
% y = img_dct2; % No need to repeat this as known DCT coefficients have
% been assigned in the last step for imgs.x01.
y(mask_neg) = lb_y(mask_neg);
BERs.x00 = sum(y(mask_neg)~=img_dct2(mask_neg))/mask_neg_num;
imgs.x00 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
% Recover the image using a naive method (randomly assigning all unknown
% sign bits to 0 or 1 with equal probability)
% y = img_dct2; % As above, no need to repeat.
randNum = randi([0 1],mask_neg_num,1);
y(mask_neg) = randNum .* ub_y(mask_neg) + (1 - randNum) .* lb_y(mask_neg);
BERs.x02 = sum(y(mask_neg)~=img_dct2(mask_neg))/mask_neg_num;
imgs.x02 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);

Ts.seconds_x0123 = toc(t1);


Ts.seconds_all = toc(t0);

disp(separator);

fprintf('\nComparison of the global sum of smoothness (i.e. the objective value):\n');
imgs_fields = {'x0' 'x00' 'x01' 'x02' 'x' 'x1' 'x2' 'x3' 'x4' 'x5' 'x6'};
strings = {'Original image' 'Recovered image (Native Method 1)' 'Recovered image (Native Method 2)' 'Recovered image (Native Method 3)' ...
           'Recovered image (CPLEX directly)' ...
           'Recovered image (Method 1)' 'Recovered image (Method 2)' 'Recovered image (Method 3)' 'Recovered image (Method 4)' ...
           'Recovered image (Method 5)' 'Recovered image (Method 6)'};
for i=1:numel(imgs_fields)
    if isfield(imgs, imgs_fields{i})
        fprintf('-%s: %d\n', strings{i}, im2smoothness(imgs.(imgs_fields{i})));
    end
end

fprintf('\nTime consumed:in seconds:\n');
fprintf('- Initial pre-processing of input arguments: %g seconds\n', Ts.seconds_prep);
fprintf('- Producing the initial condition of the model: %g seconds\n', Ts.seconds_x0123);
fprintf('- Creating the CPLEX optimization model: %g seconds\n', Ts.seconds_model);
fprintf('- Solving the CPLEX optimization model: %g seconds (%g ticks)\n', Ts.seconds_cplex, Ts.ticks_cplex);
fprintf('- Calculating the final solution(s): %g seconds\n', Ts.seconds_solutions);
fprintf('- TOTAL: %g seconds (%g%% consumed by CPLEX solver)\n', Ts.seconds_all, Ts.seconds_cplex/Ts.seconds_all*100);

fprintf('\nVisual quality of recovered images:\n');
PSNRs_all = zeros(1, 10);  % 10 = x + x00/x01/x02 + x1/x2/x3/x4/x5/x6
SSIMs_all = zeros(1, 10);
BERs_all = ones(1, 10);
PSNRs.x = psnr(imgs.x0, imgs.x);
SSIMs.x = ssim(imgs.x0, imgs.x);
PSNRs_all(1) = PSNRs.x;
SSIMs_all(1) = SSIMs.x;
BERs_all(1) = BERs.x;
fprintf('-By CPLEX directly: PSNR = %gdB, SSIM = %g\n', PSNRs.x, SSIMs.x);
x0123 = {'x00' 'x01' 'x02'};
strings = {'non-positive' 'non-negative' 'random'};
for i=1:numel(x0123)
    PSNRs.(x0123{i}) = psnr(imgs.x0, imgs.(x0123{i}));
    SSIMs.(x0123{i}) = ssim(imgs.x0, imgs.(x0123{i}));
    PSNRs_all(1+i) = PSNRs.(x0123{i});
    SSIMs_all(1+i) = SSIMs.(x0123{i});
    BERs_all(1+i) = BERs.(x0123{i});
    fprintf('-By Naive Method %d (%s): PSNR = %gdB, SSIM = %g (BER = %g%%)\n', i, strings{i}, PSNRs.(x0123{i}), SSIMs.(x0123{i}), 100*BERs.([x0123{i}]));
end
[PSNRs.x0123_max, PSNRs.x0123_maxindex] = max([PSNRs.x00 PSNRs.x01 PSNRs.x02]);
[SSIMs.x0123_max, SSIMs.x0123_maxindex] = max([SSIMs.x00 SSIMs.x01 SSIMs.x02]);
fprintf('-By the best naive method: PSNR = %gdB (Method %d), SSIM = %g (Method %d)\n', PSNRs.x0123_max, PSNRs.x0123_maxindex, SSIMs.x0123_max, SSIMs.x0123_maxindex);
x123456 = {'x1' 'x2' 'x3' 'x4', 'x5', 'x6'};
strings = {'all 0s' 'all positive' 'all negative' 'random', 'random bounding, uniform', 'random bounding, learning based'};
for i=1:numel(x123456)
    if isfield(imgs, x123456{i})
        PSNRs.(x123456{i}) = psnr(imgs.x0, imgs.(x123456{i}));
        SSIMs.(x123456{i}) = ssim(imgs.x0, imgs.(x123456{i}));
        PSNRs_all(4+i) = PSNRs.(x123456{i});
        SSIMs_all(4+i) = SSIMs.(x123456{i});
        BERs_all(4+i) = BERs.(x123456{i});
        fprintf('-By Recover Method %d (%s): PSNR = %gdB, SSIM = %g (BER = %g%%)\n', i, strings{i}, PSNRs.(x123456{i}), SSIMs.(x123456{i}), 100*BERs.([x123456{i}]));
        fprintf(' --IMPROVEMENT over the BEST naive method: \\Delta(PSNR) = %gdB, \\Delta(SSIM) = %g\n', PSNRs.(x123456{i})-PSNRs.x0123_max, SSIMs.(x123456{i})-SSIMs.x0123_max);
    end
end
x0123456 = cat(2, {'x'}, x0123, x123456);
[PSNRs.all_max, PSNRs_maxindex] = max(PSNRs_all);
PSNRs.all_maxname = x0123456{PSNRs_maxindex};
[SSIMs.all_max, SSIMs_maxindex] = max(SSIMs_all);
SSIMs.all_maxname = x0123456{SSIMs_maxindex};
[BERs.all_min, BERs_minindex] = min(BERs_all);
BERs.all_minname = x0123456{BERs_minindex};
method_names = {'CPLEX directly' 'Native Method 1' 'Native Method 2' 'Native Method 3' 'Method 1' 'Method 2' 'Method 3' 'Method 4' 'Method 5' 'Method 6'};
fprintf('-OVERALL SUMMARY:\n');
fprintf(' --The method with the best/maximum PSNR value (%g): %s\n', PSNRs.all_max, method_names{PSNRs_maxindex});
fprintf(' --The method with the best/maximum SSIM value (%g): %s\n', SSIMs.all_max, method_names{SSIMs_maxindex});
fprintf(' --The method with the best/minimum BER value (%g%%): %s\n', 100*BERs.all_min, method_names{BERs_minindex});

% modify folder_name to change the output dir
detailed_output = true; % false: output recovered images and measures
% global output_folder_name;
if ~isempty(img_filename)
    [folder_name, filename, ext] = fileparts(img_filename);
%     if ~isempty(output_folder_name)
%         folder_name = output_folder_name;
%     end
    if detailed_output
        save(fullfile(folder_name, [filename ext '_signbitLP_DCpred' num2str(DCpred_mode) '_U' num2str(round(mask_neg_num/(h/H)/(w/W))) '_T' num2str(DCTabs_min) '_relaxX' num2str(relaxX) '_relaxZ' num2str(relaxZ) '.mat']), 'imgs', 'Ts', 'imgs_dct2', 'PSNRs', 'SSIMs', 'BERs', ...
            'mask', 'DCpred_mode', 'DCTabs_min', 'recovery_method', 'DC_shift_level', 'H', 'W', 'Xmax', ...
            'twosComp', 'xmin', 'xmax', 'mask_neg_Z_num', 'cplex_solution','mask_neg_GtDCTabsMin');
    else
        save(fullfile(folder_name, [filename ext '_signbitLP_DCpred' num2str(DCpred_mode) '_U' num2str(round(mask_neg_num/(h/H)/(w/W))) '_T' num2str(DCTabs_min) '_relaxX' num2str(relaxX) '_relaxZ' num2str(relaxZ) '.mat']), 'imgs', 'Ts', 'PSNRs', 'SSIMs', 'BERs', ...
            'DCpred_mode', 'DCTabs_min', 'recovery_method', 'DC_shift_level', 'H', 'W', 'Xmax', 'twosComp', 'xmin', 'xmax', 'mask_neg_Z_num');
    end
end

% Show the figures only when there are no output arguments.
if nargout==0 && false
    if (mask_neg_Z_num>0)
        plot_images([3 4], imgs.x0, 'Original', imgs.x, 'Recovered (CPLEX)', imgs.x1, 'Recovered (Method 1)', ...
            imgs.x2, 'Recovered (Method 2)', imgs.x3, 'Recovered (Method 3)', imgs.x4, 'Recovered (Method 4)', imgs.x5, 'Recovered (Method 5)', ...
            imgs.x01, 'Recovered (Naive Method: non-negative)', imgs.x00, 'Recovered (Naive Method: non-positive)', imgs.x02, 'Recovered (Naive Method: random)');
    else % In this case no need to show imgs.x2, imgs.x3 and imgs.x4 because they are identical with imgs.x1.
        plot_images([2 4], imgs.x0, 'Original', imgs.x, 'CPLEX direct', imgs.x1, 'Recovered (Methods 1-4)', imgs.x5, 'Recovered (Method 5)', ...
            imgs.x01, 'Recovered (Naive Method: non-negative)', imgs.x00, 'Recovered (Naive Method: non-positive)', imgs.x02, 'Recovered (Naive Method: random)');
    end
    plot_dctabs2dct(imgs_dct2.y0, imgs_dct2.y, mask);
end

% Delete the log file generated by CPLEX as this is not needed.
if exist('clone1.log','file')
    delete('clone1.log');
end

end % End of the main cplex_signbit() function


function [x, y, z] = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax)
% This internal function recover x (image) and z (normal DCT coefficients)
% from y (DCT coefficients with differentially encoded DC coefficients).
% y will also be updated and returned due to the possible adjustment on x.

z = dct2blockwise_DCrecov(y, DCpred_mode, H, W);
x = idct2blockwise(z) + DC_shift_level;
x = im_postprocess(x, Xmax);
z = dct2blockwise(double(x)-DC_shift_level);
y = dct2blockwise_DCpred(z, DCpred_mode, H, W);

end % End of y2xyz() function

function x = im_postprocess(x, Xmax)
% This internal function post-processes a given image to be center shifted
% and/or its dynamic range rescaled and/or converted to uint8.

if (max(x(:))-min(x(:))>Xmax)
    % When dynamic range is larger than Xmax, rescale it.
    x = imscale(x, Xmax, 0);
else
    % When dynamic range is not larger than Xmax, centre shift it.
    x = shift2center(x, 1, 0, Xmax);
    % If the image is an 8-bit image, convert it to uint8.
    if Xmax==255
        x = uint8(x);
    end
end

end % End of im_postprocess() function
