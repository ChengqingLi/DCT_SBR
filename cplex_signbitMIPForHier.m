function [imgs, Ts, imgs_dct2, SB] = cplex_signbitMIPForHier(img_dct2, mask, DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, quan_y, quan_table, recovery_method, bDClevelshifted, H, W, Xmax, fileLabel, DCnoRelaxMask, DCLP, relaxXBound, DCoffsets, pred_method)
% [imgs, Ts, imgs_dct2, SB] = cplex_signbitMIPForHier(img_dct2, mask, DCpred_mode, DCTabs_min, timeLimit, relatGap, twosComp, relaxZBound, quan_y, quan_table, recovery_method, bDClevelshifted, H, W, Xmax, fileLabel, DCnoRelaxMask, DCLP, relaxXBound, DCoffsets, pred_method)
% This function is the subfunction invoked by cplex_signbitHierarchyDCACMIP
% to solve the problem of recovering unknown sign bits of DCT coefficients 
% of "img_dct2" (the DCT-transformed image of each region) with a mixed
% integer programming model.
%
% Input arguments:
% img_dct2: the DCT-transformed image with unknown coefficients (Note
%           img_dct2 is not necessarily the ground truth, only its absolute
%           values are meaningful)
% mask: a logical mask array showing which coefficients are unknown
%       (default: sign bits of all DC coefficients are unknown)
%       (can take a binary image or a filename as well)
% DCpred_mode: mode of DC coefficient prediction (default: 0)
%       0: no prediction
%       1, 2, 3: see dct2blockwise_DCpred
% DCTabs_min: threshold for removing DC coefficients with unknown sign bits
%             but smaller amplitude (if amplitude <= DCTabs_min, the
%             coeffient will be set to 0 and not be recovered; default: 0)
%             NB: When DCTabs_min is too large it may create too strong
%             constraints so that a feasible solution does not exist!
% timeLimit: time limit for CPLEX MIP sover. If time limit=[], the time 
%             limit is determined by CPLEX default (default: []).
% relatGap: the CPLEX relative MIP gap tolerance, if relatGap is [], 
%             relative MIP gap tolerance is determined by CPLEX default
%             (default: [])
% twosComp: whether DC is 2's complement coded, true=2's complement coded,
%           false=not 2's complement coded (default)
% relaxZBound: if set to true, z-variables are relaxed by the quantization 
%           error when twosComp is true (default: true)
% quan_y: quantized DC shifted and predicted coefficients (used only when
%           twosComp==true)
% quan_table: quantization table for JPEG files (used only when
%           twosComp==true)
% recovery_method: the method for recovering the unknown sign bits 
%                   with non-zero absolute values < DCTabs_min
%                   0: try all
%                   1: set the SB to 1
%                   2: set the SB to 0
%                   3: set the SB randomly
%                   other number: try none, leave those coefficients as 0
%                   The results of setting those coefficients to zeros
%                   is generated as CPLEX output (default)
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
%       filename. 
% DCnoRelaxMask: A mask of size [h w] specifying which (DC) coefficients will
%               require relaxation (relaxtion is required to add tolerance 
%               for DC offset individual regions in DCpred_mode 1,2,3).
%               false=need relaxtion, true=no relaxation, default: all 
%               coefficients do not need relaxation.
% DCLP: when set to true, the sign bit of DC coefficients is a floating
%               point number in [0,1]
% relaxXBound: whether the bound of x-variable is relax to [-inf, inf], 
%               true=relax, false=not relax (default)
% DCoffsets: The offsets for DC coefficients from previous regions
% pred_method: DC prediction method, used only when DCpred_mode==3, in this
%               case, when pred_method is not empty, 0: y(i,j)=z(i,j), 1:
%               y(i,j)=z(i,j)-z(i,j-W), 2: y(i,j)=z(i,j)-z(i-H,j), 3: 
%               y(i,j)=z(i,j)-0.5*(z(i-H,j)+z(i,j-W)); when pred_method
%               (default) is empty, in DCpred_mode 3, DC coefficients will
%               be calculated as if there's no DC dependency on previous 
%               regions, so the prediction methods are determine internally
%               by the positions
%
% Output arguments:
% imgs: a cell array holding a number of images recovered using different
%       methods
%       imgs.x: the image recovered directly by CPLEX
%       imgs.x1: the image recovered from recovery method 1 (see the
%       explanation for the input argument recovery_method for more
%       details, the same for the next four fields)
%       imgs.x2: the image recovered from recovery method 2
%       imgs.x3: the image recovered from recovery method 3
% Ts: the time consumed by different part of the function
%     Ts.seconds_prep: time consumed by preprocessing input arguments in seconds
%     Ts.seconds_model: time consumed by creating the CPLEX optimization model
%     Ts.seconds_cplex: time consumed by solving the CPLEX optimization model in seconds
%     Ts.ticks_cplex: deterministic time consumed by solving the CPLEX optimization model in ticks
%     Ts.seconds_solutions: time consumed by creating the final solutions in seconds
%     Ts.seconds_all: time consumed by the whole function in seconds
% imgs_dct2: DCT coefficients of some images (with differentially encoded
%            DC coefficients if DC prediction is involved)
%            imgs_dct2.y: DCT coefficients of img.x (the one recovered
%                         directly by CPLEX)
%            The above are useful for producing the figure on accuracy of
%            sign bit estimation using the optimization model.
% SB: the recovered sign bit of DCT coefficient 
%       sign bit = 1 if DCT coefficient>0
%       sign bit = 0 if DCT coefficient<=0
%       SB.x: generated by CPLEX directly
%       SB.x1: method 1
%       SB.x2: method 2
%       SB.x3: method 3
% 
% If the input image is a JPEG image (i.e. twosComp=true), then the 
% function will handle the DCT coefficients according to the quantization 
% matrix and the two's complement properly as the model will differ.
% 
% The mask can be defined to be of the same size as img_dct2, which allows
% different blocks to have different sets of unknown coefficients.
% If mask is of size HxW, it will be extended to be an array of the same
% size as the input image.
% 
% Ruiyuan Lin, Shujun Li @ University of Surrey, UK 2010-2015
% Shujun Li: http://www.hooklee.com

t0 = tic; % Start time of the whole function.

imgs = [];
Ts = struct('seconds_prep',0,'seconds_model',0,'seconds_cplex',0,'ticks_cplex',0,'seconds_solutions',0,'seconds_all',0);
imgs_dct2 = [];
SB = [];
img_filename = [];

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

xmin = 0;
if ~exist('Xmax','var')
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
if ~exist('twosComp','var')
    twosComp = false;
end
if ~exist('relaxZBound','var')
    relaxZBound = true;
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
% The correction is done before the mask is passed as an input.
% if ~bDClevelshifted
%     if DCpred_mode == 0
%         % If DC coefficients are not level shifted, then the DC coefficients'
%         % sign bits are always 1 (known), so the mask should be adjusted.
%         mask(1:H:h,1:W:w) = true;
%     elseif DCpred_mode == 1
%         mask(1:H:h,1) = true;
%     elseif DCpred_mode == 2 || DCpred_mode == 3
%         mask(1,1) = true;
%     end
% end
% if all(mask(:)) % All-true matrix
%     disp('No any sign bits are unknown, so nothing needs recovering!');
%     return;
% end

if ~exist('recovery_method','var')
    recovery_method = 4;
end
if ~exist('DCpred_mode','var')
    DCpred_mode = 0;
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
if ~exist('DCnoRelaxMask','var')
    DCnoRelaxMask = true(h,w);
end
if ~exist('DCLP','var')
    DCLP = false;
end
if ~exist('relaxXBound','var')
    relaxXBound = false;
end
if ~exist('DCoffsets','var')
    DCoffsets = zeros(h,w);
end
if ~exist('pred_method','var')
    pred_method = [];
end

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
% -- SB(1)...SB(wh) - sign bits of DCT coefficients (without DC prediction)
% -- QZ(1)...QZ(wh) - (quantized z value).*(quantization table)
% -- z(1)...z(wh) - DCT coefficients (with DC prediction)
% -- h(x(i),x(j)) - absolute values between adjacent pixels, where x(i)
%                   and x(j) are adjacent pixel values
% Objective:
% -- minimize \sum_{x(i),x(j)} h(x(i),x(j))
% Constraints:
% -- y(i)=z(i) or y(i)=z(i)-z(j) or y(i)=z(i)-(z(j1)+z(j2))/2 for different
%    DC prediction modes
% -- y(u,v) = SB(u,v)*ub_y(u,v)+(1-SB(u,v))*lb_y(u,v)
% -- x(k,l) = \sum_u\sum_v A(u,v,k,l)*z(u,v)
% -- x(i)-x(j) <= h(x(i),x(j))
% -- x(j)-x(i) <= h(x(i),x(j))
% -- y(u,v) = y*(u,v) for known y(u,v) whose value is y*(u,v)
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
% Set the display function
cplex.DisplayFunc = @redirect;

fprintf('Setting the objective...');
h_number = wh2 - w - h; % number of h-variables = adjacent pixel pairs
var_number = wh4 + h_number;
% The variables: [{x(i)}_i {SB(i)}_i {QZ(i)}_i {z(i)}_i {h(x(i),x(j)}_{x(i),x(j)}]
% Objective is 0*{x(i)}_i + 0*{SB(i)}_i + 0*{QZ(i)}_i + 0*{z(i)}_i + {h(x(i),x(j)}_{x(i),x(j)}
cplex.Model.obj = [zeros(1,wh4) ones(1,h_number)];
disp('done.');

fprintf('Setting the lower and upper bounds of variables...');
% Set lb and ub for x-, z- and h-variables.
% Note that h-variables' lb-s are tight (0) so eps is not used.
% QZ-,z-variables will be uniquely determined by SB-variables, so they 
% don't need to have lb and ub set up.
if relaxXBound
    lb = [-inf(1,wh) zeros(1,wh) -inf(1,wh2) zeros(1,h_number)];
    ub = [inf(1,wh) ones(1,wh) inf(1,wh2) inf(1,h_number)];
else
    lb = [(xmin-eps)*ones(1,wh) zeros(1,wh) -inf(1,wh2) zeros(1,h_number)];
    ub = [(xmax+eps)*ones(1,wh) ones(1,wh) inf(1,wh2) (Xmax+eps)*ones(1,h_number)];
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
        % lb and ub for SB
        lb(var_i) = (sign(img_dct2(i))>0); 
        ub(var_i) = (sign(img_dct2(i))>0);
    else % DCT coefficients with unknown signs
        if twosComp && (mod(mod(i,h),H)==1 && mod(ceil(i/h),W)==1) % JPEG image, two's complement mode.
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
                    ub_y(i) = quan_y(i)*quan_table(i); % Note as img_dct2 is not necessarily ground truth, img_dct2(i)~=quan_y(i)*quan_table(i)
                else % quan_y(i)<0
                    % Get the value of quan_y(i) without the sign bit.
                    exceptSignBit = quan_y(i) + full_range - 1;
                    % Get the lb and ub of y-variables (for future use).
                    lb_y(i) = quan_y(i)*quan_table(i);
                    ub_y(i) = (exceptSignBit+half_range)*quan_table(i);
                end
                % No need to set the lb and ub of SB
                % lb(var_i)=0 and ub(var_i)=1 is default
            end
        else
            lb_y(i) = -img_dct2_abs(i);
            ub_y(i) = img_dct2_abs(i);
            if img_dct2_abs(i)==0 %sign bit set to 0 if abs==0
                lb(var_i)=0;
                ub(var_i)=0;
            end
        end
    end
end
% the upper and lower bound for formulating the constraint matrix
% not for recovery
lb_y_con = lb_y;
ub_y_con = ub_y;
% Set the lb and ub of small y-variables to 0.
% Note that the large value between the lower bound's absolute
% value and the upper bound is used to compare with the threshold.
% The two bounds can have different absolute values when two's
% complement is used to encode differential DC coefficients.
mask_neg_ltDCTabsMin = (max(-lb_y,ub_y)<DCTabs_min) & (~(mask));
lb_y_con(mask_neg_ltDCTabsMin) = 0;
ub_y_con(mask_neg_ltDCTabsMin) = 0;
mask_neg_ltDCTabsMin_index = find(mask_neg_ltDCTabsMin);
lb(wh+mask_neg_ltDCTabsMin_index) = 0;
ub(wh+mask_neg_ltDCTabsMin_index) = 0;

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
% The number of non-zero elements in the constraints between QZ and z.
NZ_number3 = wh2;
% The maximal number of non-zero elements in the constraint matrix A.
NZ_number = NZ_number0 + NZ_number1 + NZ_number2 + NZ_number3;
% A vector of y-coordinates of all non-zero elements.
A_index_i = zeros(NZ_number,1);
% A vector of x-coordinates of all non-zero elements.
A_index_j = zeros(NZ_number,1);
% A vector of all non-zero elements.
A_value_ij = zeros(NZ_number,1);
if twosComp && relaxZBound
    cplex.Model.lhs = [zeros(wh,1); -lb_y_con(:)-DCoffsets(:); quan_table(:)/2; -inf(2*h_number,1)];
    cplex.Model.rhs = [zeros(wh,1); -lb_y_con(:)-DCoffsets(:); quan_table(:)/2; zeros(2*h_number,1)];
else
    cplex.Model.lhs = [zeros(wh,1); -lb_y_con(:)-DCoffsets(:); zeros(wh,1); -inf(2*h_number,1)];
    cplex.Model.rhs = [zeros(wh,1); -lb_y_con(:)-DCoffsets(:); zeros(wh,1); zeros(2*h_number,1)];
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

% Add the constraints about DC prediction between SB- and z-variables.
% SB(i)*ub_y_con(i)+(1-SB(i))*lb_y_con(i)=QZ(i) or
% SB(i)*ub_y_con(i)+(1-SB(i))*lb_y_con(i) =QZ(i)-QZ(j) or
% SB(i)*ub_y_con(i)+(1-SB(i))*lb_y_con(i)=QZ(i)-(QZ(j1)+QZ(j2))/2
% They are actually represented here as:
% -lb_y_con(i) <= SB(i)*(ub_y_con(i)-lb_y_con(i)) - QZ(i) <= -lb_y_con(i) (DCpred_mode = 0) or
% -lb_y_con(i) <= SB(i)*(ub_y_con(i)-lb_y_con(i)) - QZ(i) + QZ(j) <= -lb_y_con(i) (DCpred_mode = 1 and 2) or
% -lb_y_con(i) <= SB(i)*(ub_y_con(i)-lb_y_con(i)) - QZ(i) + 0.5*QZ(j1)+ 0.5*QZ(j2) <= -lb_y_con(i) (DCpred_mode = 3)
fprintf('Setting the constraints about differentially encoded DC coefficients...');
Wh = W * h;
for j=1:w % Column by column
    for i=1:h % Row by row
        if DCnoRelaxMask(i,j)
            % Add the non-zero element corresponding to SB(i).
            NZ_index = NZ_index + 1;
            con_ind = con_ind + 1;
            A_index_i(NZ_index) = con_ind;
            A_index_j(NZ_index) = con_ind; % = wh + sub2ind([h w],i,j)
            A_value_ij(NZ_index) = ub_y_con(con_ind-wh)-lb_y_con(con_ind-wh);
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
                            if isempty(pred_method) || pred_method(i,j)==1
                                NZ_index = NZ_index + 1;
                                A_index_i(NZ_index) = con_ind;
                                A_index_j(NZ_index) = wh_plus_con_ind - Wh; % z(i,j-W)
                                A_value_ij(NZ_index) = 1;
                            elseif pred_method(i,j)==3
                                NZ_index = NZ_index + 1;
                                A_index_i(NZ_index) = con_ind;
                                A_index_j(NZ_index) = wh_plus_con_ind - Wh; % z(i,j-W)
                                A_value_ij(NZ_index) = 0.5;
                            end
                        elseif (j==1 && i>1) % >=2nd DC coefficients of the first block-column
                            if isempty(pred_method) || pred_method(i,j)==2
                                NZ_index = NZ_index + 1;
                                A_index_i(NZ_index) = con_ind;
                                A_index_j(NZ_index) = wh_plus_con_ind - H; % z(i-H,j)
                                A_value_ij(NZ_index) = 1;
                            elseif pred_method(i,j)==3
                                NZ_index = NZ_index + 1;
                                A_index_i(NZ_index) = con_ind;
                                A_index_j(NZ_index) = wh_plus_con_ind - H; % z(i-H,j)
                                A_value_ij(NZ_index) = 0.5;
                            end
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
        else
            con_ind = con_ind + 1;
            cplex.Model.lhs(con_ind) = 0;
            cplex.Model.rhs(con_ind) = 0;
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

% Set the variable type
cplex.Model.ctype = [repmat('C',[1 wh]) repmat('B',[1 wh]) repmat('C',[1 (wh2+h_number)])];
if DCLP
    for j=1:W:w
        for i=1:H:h
            if ~mask(i,j)
                cplex.Model.ctype(wh+(j-1)*h+i) = 'C';
            end
        end
    end
end

Ts.seconds_model = toc(t1);

% Solve the optimization problem
disp('Solving the optimization problem...');
separator = repmat('=',1,80);
disp(separator);

% Get the mask for DCT coefficients with unknown signs.
mask_neg = ~mask;
mask_neg_num = sum(mask_neg(:));

% Open the log file in append mode
% if ~isempty (img_filename)
%     [folder_name, filename, ext] = fileparts(img_filename);
%     [fid,message] = fopen(fullfile(folder_name,[filename ext '_DCpred' num2str(DCpred_mode) '_U' num2str(round(mask_neg_num/(h/H)/(w/W))) '_T' num2str(DCTabs_min) '_timeLimit' num2str(timeLimit) '_relatGap' num2str(relatGap) '_relaxXBound' num2str(relaxXBound) '_relaxZBound' num2str(relaxZBound) '_timeLimit' num2str(timeLimit) '_signbitMIPcplex.log']), 'a');
%     if fid < 0
%       disp(message)
%     end
% elseif exist ('fileLabel', 'var')
%     [folder_name, filename, ext] = fileparts(fileLabel);
%     [fid,message] = fopen(fullfile(folder_name,[filename ext '_DCpred' num2str(DCpred_mode) '_U' num2str(round(mask_neg_num/(h/H)/(w/W))) '_T' num2str(DCTabs_min) '_timeLimit' num2str(timeLimit) '_relatGap' num2str(relatGap) '_relaxXBound' num2str(relaxXBound) '_relaxZBound' num2str(relaxZBound) '_timeLimit' num2str(timeLimit) '_signbitMIPcplex.log']), 'a');
%     if fid < 0
%       disp(message)
%     end
% end
% remember to fclose(fid);

% Set the display level of log files
cplex.Param.mip.display.Cur = 4;
if ~isempty(relatGap)
    cplex.Param.mip.tolerances.mipgap.Cur = relatGap;
end
if ~isempty(timeLimit)
    cplex.Param.timelimit.Cur = timeLimit;
end
cplex.solve();

% Close the file
% fclose(fid);

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

disp('Calculating the final solution(s)...');
% Get the directly recovered image by CPLEX.
imgs.x = imgs.x + DC_shift_level;
% imgs.x = im_postprocess(imgs.x, Xmax);
% Recovered DCT coefficients (with differentially encoded DC coefficients)
SB.x = cplex_solution.x(wh+1:wh2);
imgs_dct2.y = SB.x .* ub_y_con(:) + (1 - SB.x) .* lb_y_con(:);
imgs_dct2.y = reshape(imgs_dct2.y, [h w]);

% Mask of unknow DC coefficients with non-zero absolute value < DCTabs_min
mask_neg_ltDCTabsMin_NZ = mask_neg_ltDCTabsMin & (ub_y~=0);
mask_neg_ltDCTabsMin_NZ_num = sum(sum(mask_neg_ltDCTabsMin_NZ));
fprintf('The number of unknown sign bits that cannot be directly estimated: %d (%g%%).\n', mask_neg_ltDCTabsMin_NZ_num, mask_neg_ltDCTabsMin_NZ_num/mask_neg_num*100);
% The variable y will be used repeatedly for all recovery methods.
y = imgs_dct2.y;
% Method 1: For DCT coefficients whose sign bits are still unknown
% (unknow DC coefficients with non-zero absolute value < DCTabs_min),
% set their sign bits to 1 (positive).
% If mask_neg_ltDCTabsMin_NZ is empty, then imgs.x1 will be identical with imgs.x, so
% no need to do this. The same for imgs.x2
if (recovery_method==0 || recovery_method==1)
    if (mask_neg_ltDCTabsMin_NZ_num>0)
        SB.x1=SB.x;
        SB.x1(mask_neg_ltDCTabsMin_NZ) = 1;
        y(mask_neg_ltDCTabsMin_NZ) = ub_y(mask_neg_ltDCTabsMin_NZ);
        imgs.x1 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
    else
        disp('All unknown sign bits can be directly estimated, so skip Recovery Method 1!');
    end
end
% Method 2: For DCT coefficients whose sign bits are still unknown
% (unknow DC coefficients with non-zero absolute value < DCTabs_min),
% set their sign bits to 0 (non-positive).
if (recovery_method==0 || recovery_method==2)
    if (mask_neg_ltDCTabsMin_NZ_num>0)
        SB.x2=SB.x;
        SB.x2(mask_neg_ltDCTabsMin_NZ) = 0;
        y(mask_neg_ltDCTabsMin_NZ) = lb_y(mask_neg_ltDCTabsMin_NZ);
        imgs.x2 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
    else
        disp('All unknown sign bits can be directly estimated, so skip Recovery Method 2!');
    end
end
% Method 3: For DCT coefficients whose sign bits are still unknown
% (unknow DC coefficients with non-zero absolute value < DCTabs_min),
% randomly set their sign bits.
if (recovery_method==0 || recovery_method==3)
    if (mask_neg_ltDCTabsMin_NZ_num>0)
        SB.x3=SB.x;
        SB.x3(mask_neg_ltDCTabsMin_NZ) = randi([0 1],mask_neg_ltDCTabsMin_NZ_num,1);
        y(mask_neg_ltDCTabsMin_NZ) = SB.x3(mask_neg_ltDCTabsMin_NZ).*ub_y(mask_neg_ltDCTabsMin_NZ)+(1-SB.x3(mask_neg_ltDCTabsMin_NZ)).*lb_y(mask_neg_ltDCTabsMin_NZ);
        imgs.x3 = y2xyz(y, DCpred_mode, H, W, DC_shift_level, Xmax);
    else
        disp('All unknown sign bits can be directly estimated, so skip Recovery Method 3!');
    end
end

disp('done.');

Ts.seconds_solutions = toc(t1);

Ts.seconds_all = toc(t0);

disp(separator);

% if ~isempty(img_filename) || exist('fileLabel','var')
%     save(fullfile(folder_name, [filename ext '_signbitMIP_DCpred' num2str(DCpred_mode) '_U' num2str(round(mask_neg_num/(h/H)/(w/W))) '_T' num2str(DCTabs_min) '_timeLimit' num2str(timeLimit)  '_relatGap' num2str(relatGap) '_relaxXBound' num2str(relaxXBound) '_relaxZBound' num2str(relaxZBound) '.mat']), 'imgs', 'Ts', 'imgs_dct2', ...
%         'mask', 'DCpred_mode', 'DCTabs_min', 'recovery_method', 'DC_shift_level', 'H', 'W', 'Xmax', 'timeLimit', 'DCnoRelaxMask', 'DCLP', 'relaxXBound', 'relaxZBound', 'DCoffsets', ...
%         'xmin', 'xmax', 'mask_neg_ltDCTabsMin_NZ_num', 'cplex_solution', 'SB', 'twosComp');
% end
% 
% Delete the log file generated by CPLEX as this is not needed.
% if exist('clone1.log','file')
%     delete('clone1.log');
% end
% 
    function redirect(l)
    % This internal function redirect the CPLEX log output to a file 

    % Write the line of log output
%     fprintf(fid, '%s\r\n', l);

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
