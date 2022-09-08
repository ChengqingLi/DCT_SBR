function [imgs, Ts, imgs_dct2, SB, PSNRs, SSIMs, BERs] = cplex_signbitHierarchyDCACMIP(img_dct2, mask, regionSize_H, regionSize_W, DCpred_mode, DCdepend, relaxXBound, relaxZBound, DCpassByRegion, DCTabs_min, timeLimit, relatGap, DCpass_method, recovery_method, relaxXBound_LP, relaxZBound_LP, bDClevelshifted, H, W, Xmax, fileLabel)
% [imgs, Ts, imgs_dct2, SB, PSNRs, SSIMs, BERs] = cplex_signbitHierarchyDCACMIP(img_dct2, mask, regionSize_H, regionSize_W, DCpred_mode, DCdepend, relaxXBound, relaxZBound, DCpassByRegion, DCTabs_min, timeLimit, relatGap, DCpass_method, recovery_method, relaxXBound_LP, relaxZBound_LP, bDClevelshifted, H, W, Xmax, fileLabel)
% This function solves the problem of recovering unknown sign bits of DCT 
% coefficients of a DCT-transformed image "img_dct2" region by region (each
% by a mixed integer programming model).
%
% Input arguments:
% img_dct2: the DCT-transformed image with unknown coefficients (can also
%           take an image or a file name as the input which will be
%           read/converted to its blockwise 8x8 2D DCT matrix)
% mask: a logical mask array showing which coefficients are unknown
%       (default: sign bits of all DC coefficients are unknown)
%       (can take a binary image or a filename as well)
% regionSize_H, regionSize_W: height or width of each region, if the image
%       size cannot not be divided by the region height or width, region in
%       the last region-column will be of width=mod(w,regionSize_W), region
%       in the last region-row will be of height=mod(h,regionSize_H).
%       (Default: regionSize_H=64, regionSize_W=regionSize_H). Note if 
%       DCpred_mode==2 && (DCdepend==1 || DCdepend==2), regionSize_W will
%       be forced to be equal to image width
% DCpred_mode: mode of DC coefficient prediction (default: 0)
%       0: no prediction
%       1, 2, 3: see dct2blockwise_DCpred
% DCdepend: the DC dependency mode used in 1st pass to deal with the
%       missing information from the DC coefficients of previous regions
%       DCdepend=0: all DC coefficients that has dependency on previous
%       regions will be relaxed to [-inf inf], the dependency within the
%       region is maintained though.
%       DCdepend=1: as the regions are solved row by row, the region is
%       extended to be a rectangle from the upper left corner to the lower
%       right corner of the current unsolved region, the sign bits from 
%       previous solved regions are marked as known and use the previously
%       solved values as their value, i.e. the DC prediction will be 
%       maintained and all the pixel differences in this extended region 
%       will be concerned.
%       DCdepend=2: the current solved region will use the DC coefficients
%       from previously solved regions as an offset to the current DC
%       coefficent, but only the pixel differences in the current unsolved
%       region is concerned. If the input is JPEG file, as the quantization
%       error accumulates along the DC prediction path, the accumulated
%       quantization error is added to the relaxation of z-varibles when
%       relaxZBound is set to true
% relaxXBound: when set to true, the x-variables are relaxed to [-inf inf]
%       in 1st pass when DCdepend=1,2 (x-variables bounds are NOT relaxed 
%       when DCdepend=0 regardless of the value of relaxXBound.Note 
%       x-variables are relaxed to [-inf inf] in 2nd pass in all cases.
% relaxZBound: when set to true and when the input is a JPEG file, the
%       bounds of z-variables will be relaxed by the quantization errors
% DCpassByRegion: mode of LP 2nd pass. 
%       0: try all methods (default)
%       1: the brightness are aligned in a region-based manner
%       2: all blocks participate in the LP 2nd pass individually
% DCTabs_min: threshold for removing DC coefficients with unknown sign bits
%             but smaller amplitude (if amplitude <= DCTabs_min, the
%             coeffient will be set to 0 and not be recovered; default: 0)
%             NB: When DCTabs_min is too large it may create too strong
%             constraints so that a feasible solution does not exist!
% timeLimit: time limit for CPLEX MIP sover. If time limit=[], the time 
%             limit is determined by CPLEX default (default: []). Note the 
%             LP 2nd pass is also limited by this parameter, as it's passed
%             to CPLEX as a MIP model with no value-changeable binaries.
% relatGap: the CPLEX relative MIP gap tolerance, if relatGap is [], 
%             relative MIP gap tolerance is determined by CPLEX default
%             (default: []). Note the LP 2nd pass is not limited by this 
%              parameter
% DCpass_method: method of 2nd pass
%             0: try all (default)
%             1: MIP
%             2: LP
% recovery_method: the method in each MIP pass for regions or MIP/LP 2nd 
%                   pass for recovering the unknown sign bits with non-zero
%                   absolute values < DCTabs_min. Note the results will be
%                   available in the mat file saved from each pass, but
%                   only the one generated directly as CPLEX output is used
%                   in this function
%                   0: try all
%                   1: set the SB to 1
%                   2: set the SB to 0
%                   3: set the SB randomly
%                   other number: try none, leave those coefficients as 0
%                   The results of setting those coefficients to zeros
%                   is generated as CPLEX output (default)
% relaxXBound_LP: relaxX used in the pure LP method (default: true)
% relaxZBound_LP: relaxZ used in the pure LP method (default: true)
% bDClevelshifted: true = level shift DC coefficients by -round(xmax/2) as
%                  defined in JPEG standard (default)
%                  false = do not level shift DC coefficients (in this case
%                  DC coefficients' sign bits are always 1 so effectively
%                  known in all cases)
% H, W: height and width of the 2-D DCT block (default: H = 8, W = H)
% Xmax: maximal pixel value (default: 255)
%       The minimum pixel value Xmin is assumed to be always 0.
% fileLabel: the (imaginary) file path for img_dct2 to be used for naming
%       the mat data file generated. Used only when the img_dct2 is not a 
%       filename
%
% Output arguments:
% imgs: a cell array holding a number of images recovered using different
%       methods
%       imgs.x0: the original image as input
%       imgs.x: 1st pass
%       imgs.x1: 2nd MIP pass
%       imgs.x2: 2nd LP pass by blocks
%       imgs.x3: 2nd LP pass by regions
%       imgs.x00: the image recovered from a naive method (all unknown
%                 signs are set to 0)
%       imgs.x01: the image recovered from a naive method (all unknown
%                 signs are set to 1)
%       imgs.x02: the image recovered from a naive method (all unknown
%                 signs are randomly assigned to 0 or 1)
%       imgs.{'x000' 'x001' 'x002' 'x003' 'x004', 'x005', 'x006'}: pure LP
%                 methods, please refer to cplex_signbit for details
% Ts: the time consumed by different part of the function
%     Ts.seconds_prep: time consumed by preprocessing input arguments in
%     seconds (sumed up in all passes)
%     Ts.seconds_x0123: time consumed by creating initial condition(s)
%     Ts.seconds_model: time consumed by creating the Gurobi optimization
%     model (sumed up in all passes)
%     Ts.seconds_cplex: time consumed by solving the Gurobi optimization
%     model in seconds (sumed up in all passes)
%     Ts.ticks_cplex: deterministic time consumed by solving the CPLEX
%     optimization model in ticks (sumed up in all passes)
%     Ts.seconds_solutions: time consumed by creating the final solutions
%     in seconds (sumed up in all passes)
%     Ts.seconds_all: time consumed by the whole function in seconds
%     Ts.seconds_1stPass: time consumed in 1st pass
%     Ts.seconds_cplex_1stPass, Ts.ticks_cplex_1stPass: a matrix containing 
%     cplex seconds, ticks for each region
%     Ts.ticks_cplex_2ndPassMIP, Ts.seconds_cplex_2ndPassMIP: time comsumed
%     by CPLEX in 2nd MIP pass
%     Ts.seconds_2ndPassMIP: total time in 2nd MIP pass
%     Ts.seconds_cplex_2ndPassLP, Ts.seconds_ticks_2ndPassLP: time comsumed
%     by CPLEX in 2nd LP pass by blocks
%     Ts.seconds_cplex_2ndPassLP_region, Ts.seconds_ticks_2ndPassLP_region:
%     time consumed by CPLEX in 2nd LP pass by regions 
%     Ts.seconds_2ndPassLP: total time in 2nd LP pass
%     Ts.seconds_cplex_LP, Ts.ticks_cplex_LP: time consumed by CPLEX in
%     pure LP methods
% imgs_dct2: DCT coefficients of some images (with differentially encoded
%            DC coefficients if DC prediction is involved)
%            imgs_dct2.y0: DCT coefficients of the input (original) image
%            imgs_dct2.y: 1st pass
%            imgs_dct2.y1: MIP 2nd pass
%            imgs_dct2.y2: LP 2nd pass by blocks
%            The above are useful for producing the figure on accuracy of
%            sign bit estimation using the optimization model.
% SB: the recovered sign bit of DCT coefficient 
%       sign bit = 1 if DCT coefficient>0
%       sign bit = 0 if DCT coefficient<=0
%       SB.x: 1st pass
%       SB.x1: MIP 2nd pass
%       SB.x2: LP 2nd pass by blocks(SB is a floating point number with 0<=SB<=1 in 
%       this case)
% PSNRs, SSIMs: indicators of visual quality of recoverd images as measured
%               by PSNR and SSIM
% BERs: Bit Error Rates of different recovery methods
%       BER = Number of wrong bits / Total number of bits
%       Note BERs are calculated with SB even for LP 2nd passes
%
% If the input image is a JPEG image (when a filename is given as the first
% argument), then the function will handle the DCT coefficients according
% to the quantization matrix and the two's complement properly as the model
% will differ.
% 
% The mask can be defined to be of the same size as img_dct2, which allows
% different blocks to have different sets of unknown coefficients.
% If mask is of size HxW, it will be extended to be an array of the same
% size as the input image.
%
% Example: 
%     Region size 64x64, DCpred_mode=1, DCdepend=1, x bounds relaxed for 1st pass, try both DC LP pass by blocks and by regions, DCTabs_min=0, time limit=10min (use CPLEX default relative gap tolerance):
%     [imgs, Ts, imgs_dct2, SB, PSNRs, SSIMs, BERs] = cplex_signbitHierarchyDCACMIP('images\cameraman.pgm', 'masks\mask_9.pgm', 64, 64, 1, 1, true, false, 0, 0, 600);
%     Note: relaxZBound is not relevant to non-jpeg image, its can be set whatever values in that case
%
% Ruiyuan Lin, Shujun Li @ University of Surrey, UK 2010-2015
% Shujun Li: http://www.hooklee.com
t0 = tic; % Start time of the whole function.

imgs = [];
Ts = struct('seconds_prep',0,'seconds_model',0,'seconds_cplex',0,'ticks_cplex',0, ...
        'seconds_solutions',0,'seconds_all',0,'seconds_x0123',0,'seconds_1stPass',0, ...
        'ticks_cplex_1stPass',[],'seconds_cplex_1stPass',[],'ticks_cplex_2ndPassMIP',0, ...
        'seconds_cplex_2ndPassMIP',0,'seconds_2ndPassMIP',0,'seconds_cplex_2ndPassLP',0, ...
        'seconds_ticks_2ndPassLP',0,'seconds_2ndPassLP',0,'seconds_cplex_LP',0,'ticks_cplex_LP',0, ...
        'seconds_cplex_2ndPassLP_region',0,'seconds_ticks_2ndPassLP_region',0);
imgs_dct2 = [];

SB = [];
PSNRs = [];
SSIMs = [];
BERs = [];

img_filename = [];

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

if ~exist('Xmax','var')
    Xmax = 255;
end
if ~exist('bDClevelshifted','var')
    bDClevelshifted = true;
end
if bDClevelshifted
    DC_shift_level = round(Xmax/2); % e.g. 255 => 128 (JPEG standard)
else
    DC_shift_level = 0;
end
if ~exist('H','var')
    H = 8;
end
if ~exist('W','var')
    W = H;
end
if ~exist('regionSize_H','var')
    regionSize_H = 64;
end
if ~exist('regionSize_W','var')
    regionSize_W = regionSize_H;
end
if ~exist('DCpass_method','var')
    DCpass_method = 0; 
end
if ~exist('relaxXBound_LP','var')
    relaxXBound_LP = true;
end
if ~exist('relaxZBound_LP','var')
    relaxZBound_LP = true;
end
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
if (w<regionSize_W || h<regionSize_H)
    disp('The region size is larger than the image size!');
    return;
end
if (mod(regionSize_W,W)>0 || mod(regionSize_H,H)>0)
    disp('The region size does not match the block size!');
    return;
end
if DCpred_mode==2 && (DCdepend==1 || DCdepend==2)
    if regionSize_W~=w
        disp('DCpred_mode 2, DCdpend 1 or 2, forcing region width equals image width!');
        regionSize_W = w;
    end
end
mask = process_mask(mask, h, w);
if ~bDClevelshifted
    if DCpred_mode==0
        % If DC coefficients are not level shifted, then the DC coefficients'
        % sign bits are always 1 (known), so the mask should be adjusted.
        mask(1:H:h,1:W:w) = true;
    elseif DCpred_mode==1
        mask(1:H:h,1) = true;
    elseif DCpred_mode==2 || DCpred_mode==3
        mask(1,1) = true;
    end
end
if all(mask(:)) % All-true matrix
    disp('No any sign bits are unknown, so nothing needs recovering!');
    return;
end
mask_neg = ~mask;
mask_neg_num = sum(mask_neg(:));
U = round(mask_neg_num/(h/H)/(w/W)); % used for naming the mat files

if ~exist('recovery_method','var')
   recovery_method = 4;
end
if ~exist('DCpred_mode','var')
   DCpred_mode = 0;
end
if ~exist('DCdepend','var')
    DCdepend = 0;
end
if ~exist('relaxXBound','var')
    relaxXBound = false;
end
if ~exist('relaxZBound','var')
    relaxZBound = false;
end
if ~exist('DCpassByRegion','var')
    DCpassByRegion = 0;
end
if ~exist('DCTabs_min','var')
   DCTabs_min = 0; 
end
if ~exist('timeLimit','var')
   timeLimit = []; 
end
if ~exist('relatGap','var')
   relatGap = [];
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
    imgs.x0 = uint8(imgs.x0);
end
imgs_dct2.y0 = img_dct2;
SB_GT = (sign(img_dct2)>0); 

seg_num_h = ceil(h/regionSize_H);
seg_num_w = ceil(w/regionSize_W);
segmentFileNames = cell(seg_num_h,seg_num_w);
imgsSeg = cell(seg_num_h,seg_num_w);
TsSeg = cell(seg_num_h,seg_num_w);
imgs_dct2Seg = cell(seg_num_h,seg_num_w);
SBSeg = cell(seg_num_h,seg_num_w);
Ts.seconds_cplex_1stPass = zeros(seg_num_h,seg_num_w);
Ts.ticks_cplex_1stPass = zeros(seg_num_h,seg_num_w);
        
if ~isempty (img_filename)
    [folder_name, filename, ext] = fileparts(img_filename);
elseif exist ('fileLabel', 'var')
    [folder_name, filename, ext] = fileparts(fileLabel);
end

DCNoRelaxMask = true(h,w);
% if DCpred_mode == 1 || DCpred_mode == 2 || DCpred_mode == 3
%     DCNoRelaxMask(1:H:h,1:W:w) = false;
% end
% relax DC coefficients that has dependency on previous regions to [-inf, inf]
if DCdepend==0
    if DCpred_mode == 1
        DCNoRelaxMask(1:H:h,(1+regionSize_W):regionSize_W:w) = false; % the first of each row-block in each region starting from the second region-column
    elseif DCpred_mode == 2
        DCNoRelaxMask(1:H:h,1:regionSize_W:w) = false; 
        DCNoRelaxMask(1,1) = true; % the first of each row-block in each region starting except the very first DC coefficient
        if regionSize_W==w % in this case all except the first all every region have correct dependency
            DCNoRelaxMask(1:H:h,1) = true;
            DCNoRelaxMask((1+regionSize_H):regionSize_H:h,1) = false;
        end
    elseif DCpred_mode == 3
        DCNoRelaxMask(1:regionSize_H:h,1:W:w) = false; % the first row-block of each region
        DCNoRelaxMask(1:H:h,1:regionSize_W:w) = false; % the first column-block of each region
        DCNoRelaxMask(1,1:W:w) = true; 
        DCNoRelaxMask(1,(1+regionSize_W):regionSize_W:w) = false; % first row except the first of each block
        DCNoRelaxMask(1:H:h,1) = true; 
        DCNoRelaxMask((1+regionSize_H):regionSize_H:h,1) = false; % first column except the first of each block
    end
end
imgs.x = zeros(h,w); % direct combination of segments
imgs_dct2.y = img_dct2; % direct combination of segments
SB.x = zeros(h,w);
imgs_dct2.ub_y = zeros(h,w);
imgs_dct2.lb_y = zeros(h,w);
DCoffsets = zeros(h,w);
if twosComp
    DCprevQuanError = zeros(h,w);
end
if DCpred_mode==3 && DCdepend==2
    % predict current DC from the left and/or upper blocks
    pred_method = zeros(h,w);
    pred_method((1+H):H:h,(1+W):W:w) = 3;
    pred_method(1,(1+W):W:w) = 1;
    pred_method((1+H):H:h,1) = 2;
end
Ts.seconds_prep = toc(t0);
t1 = tic;
i_i = 0;

for i = 1:regionSize_H:h
    i_i = i_i + 1;
    j_i = 0;
    for j = 1:regionSize_W:w
        j_i = j_i + 1;
        segmentFileNames{i_i,j_i} = fullfile(folder_name, [filename ext '_regionSize' num2str(regionSize_H) 'x' num2str(regionSize_W) '_row' num2str(i_i) '_col' num2str(j_i) '_DCdpend' num2str(DCdepend) ext]);
        % when h % regionSize ~= 0
        if (i+regionSize_H-1)<=h
            regionSize_H_cur = regionSize_H;
        else
            regionSize_H_cur = h-i+1;
        end
        if (j+regionSize_W-1)<=w
            regionSize_W_cur = regionSize_W;
        else
            regionSize_W_cur = w-j+1;
        end
        if DCdepend==0
            if twosComp
                [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), mask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, quan_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), quan_table(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), false, false);
            else
                [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), mask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, [], [], recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), false, false);
            end
            imgs.x(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgsSeg{i_i,j_i}.x;
            imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.y; % the DC returned value have correct absolute value (as long as it absolution value is larger than DCTabs_min) as it's calculated from ub_y, lb_y directly
            imgs_dct2.ub_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.ub_y;
            imgs_dct2.lb_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.lb_y;
            SB.x(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = reshape(SBSeg{i_i,j_i}.x, [regionSize_H_cur regionSize_W_cur]);
        elseif DCdepend==1
            mask0 = true(i-1+regionSize_H_cur,j-1+regionSize_W_cur);
            mask0(i:(i-1+regionSize_H_cur),j:(j-1+regionSize_W_cur)) = mask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1));
            if twosComp
                [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(1:(i+regionSize_H_cur-1),1:(j+regionSize_W_cur-1)), mask0, DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, quan_y(1:(i+regionSize_H_cur-1),1:(j+regionSize_W_cur-1)), quan_table(1:(i+regionSize_H_cur-1),1:(j+regionSize_W_cur-1)), recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(1:(i+regionSize_H_cur-1),1:(j+regionSize_W_cur-1)), false, relaxXBound);
            else
                [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(1:(i+regionSize_H_cur-1),1:(j+regionSize_W_cur-1)), mask0, DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, [], [], recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(1:(i+regionSize_H_cur-1),1:(j+regionSize_W_cur-1)), false, relaxXBound);
            end
            imgs.x(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgsSeg{i_i,j_i}.x(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1));
            imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)); 
            imgs_dct2.ub_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.ub_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1));
            imgs_dct2.lb_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.lb_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1));
            SB0 = reshape(SBSeg{i_i,j_i}.x, [(i+regionSize_H_cur-1) (j+regionSize_W_cur-1)]);
            SB.x(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = SB0(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1));
        elseif DCdepend==2
            if DCpred_mode~=3
                if DCpred_mode==1
                    if j>regionSize_W
                       z_prevRegions = dct2blockwise_DCrecov(imgs_dct2.y(i:(i+regionSize_H_cur-1),1:(j-1)), DCpred_mode, H, W);
                       DCoffsets(i:H:(i+regionSize_H_cur-1),j) = z_prevRegions(1:H:regionSize_H_cur,j-W);
                    end
                elseif DCpred_mode==2
                    if i>1 % Except the first region, note the regionSize_W is forced to be equal to w in DCpred_mode2 and DCdepend2
                        z_aboveRegions = dct2blockwise_DCrecov(imgs_dct2.y(1:(i-1),1:w), DCpred_mode, H, W); 
                        DCoffsets(i,j) = z_aboveRegions(i-H,w-W+1);
                        if twosComp
                            quanError_aboveRegions = dct2blockwise_DCrecov(quan_table(1:(i-1),1:w)/2, DCpred_mode, H, W); 
                            DCprevQuanError(i,j) = quanError_aboveRegions(i-H,w-W+1);
                        end
                    end
                end
                if twosComp
                    [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), mask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, quan_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), quan_table(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), false, relaxXBound, DCoffsets(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)));
                else
                    [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), mask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, [], [], recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), false, relaxXBound, DCoffsets(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)));
                end
            else % DCpred_mode==3
                if j>regionSize_W
                   z_leftRegions = dct2blockwise_DCrecov(imgs_dct2.y(1:(i+regionSize_H_cur-1),1:(j-1)), DCpred_mode, H, W);
                   DCoffsets(i:H:(i+regionSize_H_cur-1),j) = z_leftRegions(i:H:(i+regionSize_H_cur-1),j-W)/2;
                   if twosComp
                      quanError_leftRegions = dct2blockwise_DCrecov(quan_table(1:(i+regionSize_H_cur-1),1:(j-1))/2, DCpred_mode, H, W);
                      DCprevQuanError(i:H:(i+regionSize_H_cur-1),j) = quanError_leftRegions(i:H:(i+regionSize_H_cur-1),j-W)/2;
                   end
                   if pred_method(i,j)==1
                       assert(i==1,'i~=1, pred_method wrong!')
                       DCoffsets(i,j) = z_leftRegions(i,j-W);
                       if twosComp
                          DCprevQuanError(i,j) = quanError_leftRegions(i,j-W); 
                       end
                   end
                end
                if i>regionSize_H
                   z_aboveRegions = dct2blockwise_DCrecov(imgs_dct2.y(1:(i-1),1:(j+regionSize_W_cur-1)), DCpred_mode, H, W); 
                   DCoffsets(i,j:W:(j+regionSize_W_cur-1)) = DCoffsets(i,j:W:(j+regionSize_W_cur-1)) + z_aboveRegions(i-H,j:W:(j+regionSize_W_cur-1))/2;
                   if twosComp
                      quanError_aboveRegions = dct2blockwise_DCrecov(quan_table(1:(i-1),1:(j+regionSize_W_cur-1))/2, DCpred_mode, H, W);
                      DCprevQuanError(i,j:W:(j+regionSize_W_cur-1)) = DCprevQuanError(i,j:W:(j+regionSize_W_cur-1)) + quanError_aboveRegions(i-H,j:W:(j+regionSize_W_cur-1))/2;
                   end
                   if pred_method(i,j)==2 
                       assert(j==1,'j~=1, pred_method wrong!')
                       DCoffsets(i,j) = z_aboveRegions(i-H,j);
                       if twosComp
                           DCprevQuanError(i,j) = quanError_aboveRegions(i-H,j);
                       end
                   end
                end
                if twosComp
                    [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), mask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, quan_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), quan_table(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), false, relaxXBound, DCoffsets(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), pred_method(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)));                   
                else
                    [imgsSeg{i_i,j_i}, TsSeg{i_i,j_i}, imgs_dct2Seg{i_i,j_i}, SBSeg{i_i,j_i}] = cplex_signbitMIPForHier(imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), mask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, [], [], recovery_method, bDClevelshifted, H, W, Xmax, segmentFileNames{i_i,j_i}, DCNoRelaxMask(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), false, relaxXBound, DCoffsets(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)), pred_method(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)));
                end
            end
            imgs.x(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgsSeg{i_i,j_i}.x;
            imgs_dct2.y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.y; 
            imgs_dct2.ub_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.ub_y;
            imgs_dct2.lb_y(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = imgs_dct2Seg{i_i,j_i}.lb_y;
            SB.x(i:(i+regionSize_H_cur-1),j:(j+regionSize_W_cur-1)) = reshape(SBSeg{i_i,j_i}.x, [regionSize_H_cur regionSize_W_cur]);
        end
        Ts.seconds_prep = Ts.seconds_prep + TsSeg{i_i,j_i}.seconds_prep;
        Ts.seconds_model = Ts.seconds_model + TsSeg{i_i,j_i}.seconds_model;
        Ts.seconds_cplex = Ts.seconds_cplex + TsSeg{i_i,j_i}.seconds_cplex;
        Ts.ticks_cplex = Ts.ticks_cplex + TsSeg{i_i,j_i}.ticks_cplex;
        Ts.seconds_cplex_1stPass(i_i,j_i) = TsSeg{i_i,j_i}.seconds_cplex;
        Ts.ticks_cplex_1stPass(i_i,j_i) = TsSeg{i_i,j_i}.ticks_cplex;
        Ts.seconds_solutions = Ts.seconds_solutions + TsSeg{i_i,j_i}.seconds_solutions;
    end
end
% Note the output imgs.x of each segment is not pass through im_postprocess
% inside the function so as not to create different DC offset for different
% blocks
directCombRegions_x = imgs.x;
imgs.x = im_postprocess(imgs.x, Xmax);
mask_neg_GtDCTabsMin = mask_neg & (max(-imgs_dct2.lb_y,imgs_dct2.ub_y) >= DCTabs_min);
mask_neg_GtDCTabsMin_num = sum(mask_neg_GtDCTabsMin(:));
BERs.x = sum(SB.x(mask_neg_GtDCTabsMin)~=(SB_GT(mask_neg_GtDCTabsMin)))/mask_neg_GtDCTabsMin_num;
Ts.seconds_1stPass = toc(t1);

% Further process the DC coefficients
segCombFileName = fullfile(folder_name, [filename ext '_MIPDC_regionSize' num2str(regionSize_H) 'x' num2str(regionSize_W) '_DCdpend' num2str(DCdepend) '_SegCombination' ext]);
mask_DCAC = true(h,w);
% if DCpassByRegion && DCpred_mode~=0
%     if DCpred_mode == 1
%         mask_DCAC(1:H:h,1:regionSize_W:w) = false; % the first of each row-block in each region
%     elseif DCpred_mode == 2
%         if regionSize_W==w 
%             mask_DCAC(1:regionSize_H:h,1) = false;
%         else
%             mask_DCAC(1:H:h,1:regionSize_W:w) = false; 
%         end
%     elseif DCpred_mode == 3
%         mask_DCAC(1:regionSize_H:h,1:W:w) = false; % the first row-block of each region
%         mask_DCAC(1:H:h,1:regionSize_W:w) = false; % the first column-block of each region
%         mask_DCAC(1,1:W:w) = true; 
%         mask_DCAC(1,1:regionSize_W:w) = false; % first row except the first of each block
%         mask_DCAC(1:H:h,1) = true; 
%         mask_DCAC(1:regionSize_H:h,1) = false; % first column except the first of each block
%     end
% else
%     mask_DCAC(1:H:h,1:W:w) = false;
% end
mask_DCAC(1:H:h,1:W:w) = false;
mask2 = mask | mask_DCAC;
% if all(mask2(:)) % All-true matrix
%     disp('No any DC sign bits are unknown, so skip the global DC recovery step!');
%     return;
% end
t1 = tic;
relaxXBound2 = true; % Always relax x-variable bounds for 2nd pass
if DCpass_method==0 || DCpass_method==1
    if twosComp
        [imgs_1, Ts_1, imgs_dct2_1, SB_1] = cplex_signbitMIPForHier(imgs_dct2.y, mask2, DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, quan_y, quan_table, recovery_method, bDClevelshifted, H, W, Xmax, segCombFileName, true(h,w), false, relaxXBound2);
    else
        [imgs_1, Ts_1, imgs_dct2_1, SB_1] = cplex_signbitMIPForHier(imgs_dct2.y, mask2, DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, [], [], recovery_method, bDClevelshifted, H, W, Xmax, segCombFileName, true(h,w), false, relaxXBound2);
    end
    imgs.x1 = imgs_1.x;
    imgs.x1 = im_postprocess(imgs.x1,Xmax);
    imgs_dct2.y1 = imgs_dct2_1.y;
    SB.x1 = reshape(SB_1.x, [h w]);
    Ts.seconds_prep = Ts.seconds_prep + Ts_1.seconds_prep;
    Ts.seconds_model = Ts.seconds_model + Ts_1.seconds_model;
    Ts.seconds_cplex = Ts.seconds_cplex + Ts_1.seconds_cplex;
    Ts.ticks_cplex = Ts.ticks_cplex + Ts_1.ticks_cplex;
    Ts.seconds_cplex_2ndPassMIP = Ts_1.seconds_cplex;
    Ts.ticks_cplex_2ndPassMIP = Ts_1.ticks_cplex;
    Ts.seconds_solutions = Ts.seconds_solutions + Ts_1.seconds_solutions;
    BERs.x1 = sum(SB.x1(mask_neg_GtDCTabsMin)~=(SB_GT(mask_neg_GtDCTabsMin)))/mask_neg_GtDCTabsMin_num;
end
Ts.seconds_2ndPassMIP = toc(t1);

t1 = tic;
if DCpass_method==0 || DCpass_method==2
    if DCpassByRegion==0 || DCpassByRegion==1
        [imgs_3, Ts_3] = alignBrightness(directCombRegions_x, regionSize_H, regionSize_W);
        imgs_dct2.y3 = dct2blockwise(imgs_3.x,H,W);
        SB.x3 = (imgs_dct2.y3>0);
        BERs.x3 = sum(SB.x3(mask_neg_GtDCTabsMin)~=SB_GT(mask_neg_GtDCTabsMin))/mask_neg_GtDCTabsMin_num;
        imgs.x3 = imgs_3.x;
        imgs.x3 = im_postprocess(imgs.x3,Xmax);
        Ts.seconds_prep = Ts.seconds_prep + Ts_3.seconds_prep;
        Ts.seconds_model = Ts.seconds_model + Ts_3.seconds_model;
        Ts.seconds_cplex = Ts.seconds_cplex + Ts_3.seconds_cplex;
        Ts.ticks_cplex = Ts.ticks_cplex + Ts_3.ticks_cplex;
        Ts.seconds_cplex_2ndPassLP_region = Ts_3.seconds_cplex;
        Ts.seconds_ticks_2ndPassLP_region = Ts_3.ticks_cplex;
    end
    if DCpassByRegion==0 || DCpassByRegion==2
        [imgs_2, Ts_2] = alignBrightness(directCombRegions_x, H, W);
        imgs_dct2.y2 = dct2blockwise(imgs_2.x,H,W);
        SB.x2 = (imgs_dct2.y2>0);
        BERs.x2 = sum(SB.x2(mask_neg_GtDCTabsMin)~=SB_GT(mask_neg_GtDCTabsMin))/mask_neg_GtDCTabsMin_num;
        imgs.x2 = imgs_2.x;
        imgs.x2 = im_postprocess(imgs.x2,Xmax);
        Ts.seconds_prep = Ts.seconds_prep + Ts_2.seconds_prep;
        Ts.seconds_model = Ts.seconds_model + Ts_2.seconds_model;
        Ts.seconds_cplex = Ts.seconds_cplex + Ts_2.seconds_cplex;
        Ts.ticks_cplex = Ts.ticks_cplex + Ts_2.ticks_cplex;
        Ts.seconds_cplex_2ndPassLP = Ts_2.seconds_cplex;
        Ts.seconds_ticks_2ndPassLP = Ts_2.ticks_cplex;
    end
end
Ts.seconds_2ndPassLP = toc(t1);

% % LP + naive methods
% if ~isempty(img_filename) % in the case of JPEG file, the input should be filename
%     [imgs_LP, Ts_LP, imgs_dct2_LP, PSNRs_LP, SSIMs_LP, BERs_LP] = cplex_signbit(img_filename, mask, DCpred_mode, DCTabs_min, relaxXBound_LP, relaxZBound_LP, 0, bDClevelshifted, H, W, Xmax);
% else
%     [imgs_LP, Ts_LP, imgs_dct2_LP, PSNRs_LP, SSIMs_LP, BERs_LP] = cplex_signbit(img_dct2, mask, DCpred_mode, DCTabs_min, relaxXBound_LP, relaxZBound_LP, 0, bDClevelshifted, H, W, Xmax);
% end
% % LP fields
% xx123456 = {'x' 'x1' 'x2' 'x3' 'x4', 'x5', 'x6'};
x000123456 = {'x000' 'x001' 'x002' 'x003' 'x004', 'x005', 'x006'};
% for i=1:numel(xx123456)
%     if isfield(imgs_LP, xx123456{i})
%         imgs.(x000123456{i}) = imgs_LP.(xx123456{i});
%         PSNRs.(x000123456{i}) = PSNRs_LP.(xx123456{i});
%         SSIMs.(x000123456{i}) = SSIMs_LP.(xx123456{i});
%         BERs.(x000123456{i}) = BERs_LP.(xx123456{i});
%     end
% end
% imgs_dct2.y000 = imgs_dct2_LP.y;
% Ts.seconds_cplex_LP = Ts_LP.seconds_cplex;
% Ts.ticks_cplex_LP = Ts_LP.ticks_cplex;
% % naive method fields
x0123 = {'x00' 'x01' 'x02'};
% for i=1:numel(x0123)
%     imgs.(x0123{i}) = imgs_LP.(x0123{i});
%     PSNRs.(x0123{i}) = PSNRs_LP.(x0123{i});
%     SSIMs.(x0123{i}) = SSIMs_LP.(x0123{i});
%     BERs.(x0123{i}) = BERs_LP.(x0123{i});
% end
% PSNRs.x0123_max = PSNRs_LP.x0123_max;
% SSIMs.x0123_max = SSIMs_LP.x0123_max;
% Ts.seconds_x0123 = Ts_LP.seconds_x0123;
% % Comparison results
% SSIMs.LP_naive_max = SSIMs_LP.all_max;
% PSNRs.LP_naive_max = PSNRs_LP.all_max;
% BERs.LP_naive_min = BERs_LP.all_min;

Ts.seconds_all = toc(t0);

separator = repmat('=',1,80);
disp(separator);
fprintf('\nComparison of the global sum of smoothness (i.e. the objective value):\n');
imgs_fields = {'x0' 'x' 'x1' 'x2' 'x3' 'x00' 'x01' 'x02' 'x000' 'x001' 'x002' 'x003' 'x004', 'x005', 'x006'};
strings = {'Original image' 'Recovered image (direct combination of segments)' 'Recovered image (after the DC-only MIP pass)' 'Recovered image (after the DC-only LP pass)' 'Recovered image (after the DC-only LP pass by region)' 'Recovered image (Native Method 1)' 'Recovered image (Native Method 2)' 'Recovered image (Native Method 3)' 'LP:CPLEX direct' 'LP: all 0s' 'LP: all positive' 'LP: all negative' 'LP: random', 'LP: random bounding, uniform', 'LP: random bounding, learning based'};
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
PSNRs_all = zeros(1, 14);  % 13 = x +x1/x2/x3 +x00/x01/x02 +x000/x001/x002/x003/x004/x005/x006
SSIMs_all = zeros(1, 14);
BERs_all = ones(1, 14);
% for i=1:numel(x0123)
%     PSNRs_all(i) = PSNRs.(x0123{i});
%     SSIMs_all(i) = SSIMs.(x0123{i});
%     BERs_all(i) = BERs.(x0123{i});
% end
% for i=1:numel(x000123456)
%     if isfield(imgs, x000123456{i})
%         PSNRs_all(i+3) = PSNRs.(x000123456{i});
%         SSIMs_all(i+3) = SSIMs.(x000123456{i});
%         BERs_all(i+3) = BERs.(x000123456{i});
%     end
% end
xx123 = {'x' 'x1' 'x2' 'x3'};
strings = {'direct combination' 'after DC-only pass (MIP)' 'after DC only pass (LP)' 'after DC only pass (LP by region)'};
for i=1:numel(xx123)
    if isfield(imgs, xx123{i})
        PSNRs.(xx123{i}) = psnr(imgs.x0, imgs.(xx123{i}));
        SSIMs.(xx123{i}) = ssim(imgs.x0, imgs.(xx123{i}));
        PSNRs_all(i+10) = PSNRs.(xx123{i});
        SSIMs_all(i+10) = SSIMs.(xx123{i});
        BERs_all(i+10) = BERs.(xx123{i});
        fprintf('-By Hierarchical Recover Method %d (%s): PSNR = %gdB, SSIM = %g (BER = %g%%)\n', i, strings{i}, PSNRs.(xx123{i}), SSIMs.(xx123{i}), 100*BERs.([xx123{i}]));
    end
end
xx000123456 = cat(2, x0123, x000123456, xx123);
[PSNRs.all_max, PSNRs_maxindex] = max(PSNRs_all);
PSNRs.all_maxname = xx000123456{PSNRs_maxindex};
[SSIMs.all_max, SSIMs_maxindex] = max(SSIMs_all);
SSIMs.all_maxname = xx000123456{SSIMs_maxindex};
[BERs.all_min, BERs_minindex] = min(BERs_all);
BERs.all_minname = xx000123456{BERs_minindex};
method_names = {'Native Method 1' 'Native Method 2' 'Native Method 3' 'CPLEX directly' 'LP: Method 1' 'LP: Method 2' 'LP: Method 3' 'LP: Method 4' 'LP: Method 5' 'LP: Method 6' 'Hier: direct combination' 'Hier: after DC-only pass (MIP)' 'Hier: after DC only pass (LP)' 'Hier: after DC only pass (LP BY REGION)'};
fprintf('-OVERALL SUMMARY:\n');
fprintf(' --The method with the best/maximum PSNR value (%g): %s\n', PSNRs.all_max, method_names{PSNRs_maxindex});
fprintf(' --The method with the best/maximum SSIM value (%g): %s\n', SSIMs.all_max, method_names{SSIMs_maxindex});
fprintf(' --The method with the best/minimum BER value (%g%%): %s\n', 100*BERs.all_min, method_names{BERs_minindex});

% global output_folder_name;
% if ~isempty(output_folder_name)
%     folder_name = output_folder_name;
% end
if ~isempty(img_filename) || exist('fileLabel','var')
    save(fullfile(folder_name, [filename ext '_signbitHMIP_sz' num2str(regionSize_H) 'x' num2str(regionSize_W) '_P' num2str(DCpred_mode) '_D' num2str(DCdepend) '_U' num2str(U) '_T' num2str(DCTabs_min) '_tL' num2str(timeLimit) '_gap' num2str(relatGap) '_rX' num2str(relaxXBound) '_rZ' num2str(relaxZBound) '_DCpassByRegion' num2str(DCpassByRegion) '.mat']), 'imgs', 'Ts', 'imgs_dct2', 'SB', 'PSNRs', 'SSIMs', 'BERs', ...
        'mask', 'DCpred_mode', 'DCTabs_min', 'recovery_method', 'DC_shift_level', 'H', 'W', 'Xmax', 'timeLimit', ...
        'imgsSeg', 'TsSeg', 'imgs_dct2Seg', 'SBSeg', 'relaxXBound', 'relaxXBound2');
end

% Show the figures only when there are no output arguments.
if nargout==0
    if ~isempty(img_filename) || exist('fileLabel','var')
        figure('Name', [filename ext '_signbitHierarchyMIP_regionSize' num2str(regionSize_H) 'x' num2str(regionSize_W) '_DCpred' num2str(DCpred_mode) '_DCdpend' num2str(DCdepend) '_U' num2str(U) '_T' num2str(DCTabs_min) '_timeLimit' num2str(timeLimit) '_relaxXBound' num2str(relaxXBound) '_relaxZBound' num2str(relaxZBound)]);
    else
        figure;
    end
    if isfield(imgs,'x1') && isfield(imgs,'x2') && isfield(imgs,'x3')
    plot_images([2 3], imgs.x0, 'Original', imgs.x, 'Direct combination', ...
            imgs.x1, 'After the DC-only MIP pass', imgs.x2, 'After the DC-only LP pass', ...
            imgs.x3, 'After the DC-only LP pass (by region)');
    elseif isfield(imgs,'x1') && isfield(imgs,'x2')
    plot_images([2 2], imgs.x0, 'Original', imgs.x, 'Direct combination', ...
            imgs.x1, 'After the DC-only MIP pass', imgs.x2, 'After the DC-only LP pass');        
    elseif isfield(imgs,'x1') && isfield(imgs,'x3')
    plot_images([2 2], imgs.x0, 'Original', imgs.x, 'Direct combination', ...
            imgs.x1, 'After the DC-only MIP pass', imgs.x3, 'After the DC-only LP pass (by region)');
    end
%     plot_dctabs2dct(imgs_dct2.y0, imgs_dct2.y, mask);
%     plot_dctabs2dct(imgs_dct2.y0, imgs_dct2.y1, mask);
%     plot_dctabs2dct(imgs_dct2.y0, imgs_dct2.y2, mask);
end

end

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
