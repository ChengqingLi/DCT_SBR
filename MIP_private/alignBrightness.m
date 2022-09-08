function [imgs, Ts] = alignBrightness(img, regionSize_H, regionSize_W)
% [imgs, Ts] = alignBrightness(img, regionSize_H, regionSize_W)
% This function invokes IBM ILOG CPLEX's MATLAB Cplex class to solve the
% problem of aligning the brightness of each region of img
%
% Input arguments:
% img: input image
% regionSize_H, regionSize_W: height or width of each region
%
% Output arguments:
% imgs: a cell array holding a number of images recovered using different
%       methods
%       imgs.x: the image recovered directly by CPLEX (in double)
% Ts: the time consumed by different part of the function
%     Ts.seconds_prep: time consumed by preprocessing input arguments in seconds
%     Ts.seconds_model: time consumed by creating the CPLEX optimization model
%     Ts.seconds_cplex: time consumed by solving the CPLEX optimization model in seconds
%     Ts,ticks_cplex: deterministic time consumed by solving the CPLEX optimization model in ticks
%     Ts.seconds_all: time consumed by the whole function in seconds
%
% Ruiyuan Lin, Shujun Li @ University of Surrey, UK 2010-2015
% Shujun Li: http://www.hooklee.com

t0 = tic; % Start time of the whole function.

imgs = [];
Ts = struct('seconds_prep',0,'seconds_model',0,'seconds_cplex',0,'ticks_cplex',0,'seconds_all',0);

if nargin<1
    disp('At least one input argument is needed!');
    return;
end

img = double(img);

h = size(img,1);
w = size(img,2);

if mod(w,regionSize_W)>0 || mod(h,regionSize_H)>0
    disp('Input image does not match the region size!')
    return;
end

wh = w*h;
wh2 = wh*2;

Ts.seconds_prep = toc(t0);

t1 = tic;

% ==================The 2-D DCT sign bits recovery model===================
% Variables:
% -- x(1)...x(wh) - pixel values
% -- b(1)...b(region_num) - brightness adjustment of each region
% -- h(x(i),x(j)) - absolute values between adjacent pixels along boundary,
%                   where x(i) and x(j) are adjacent pixel values
% Objective:
% -- minimize \sum_{x(i),x(j)} h(x(i),x(j))
% Constraints:
% -- x(i) = img(i) + b(current region)
% -- x(i)-x(j) <= h(x(i),x(j))
% -- x(j)-x(i) <= h(x(i),x(j))
% An equality x=a is represented by an inequality a<=x<=a.
% Lower and upper bounds:
% -- -inf <= x(i) <= inf
% -- -inf <= b(i) <=inf
% -- 0 <= h(x(i),x(j)) <= inf
% =========================================================================

% Error tolerance, a value smaller than 0.5 is sufficient to make sure
% rounding will work properly for x-variables (pixel values).
eps = 0.499;

cplex = Cplex(mfilename);
% cplex.Model.sense = 'minimize'; % minimization is the default.
% cplex.Param.lpmethod.Cur = 0; % 'auto' algorithm selection
% Set the display function
% cplex.DisplayFunc = @redirect;

region_num_H = h/regionSize_H;
region_num_W = w/regionSize_W;
region_num = region_num_H * region_num_W;

fprintf('Setting the objective...');
h_number = wh2 - w - h; % number of h-variables = adjacent pixel pairs
var_number = wh + region_num + h_number;
% The variables: [{x(i)}_i {b(i)}_i {h(x(i),x(j)}_{x(i),x(j)}]
% Objective is 0*{x(i)}_i + 0*{b(i)}_i + {h(x(i),x(j)}_{x(i),x(j)}
cplex.Model.obj = [zeros(1,wh+region_num) ones(1,h_number)];
disp('done.');

fprintf('Setting the lower and upper bounds of variables...');
% Set lb and ub for x-, b- and h-variables.
% Note that h-variables' lb-s are tight (0) so eps is not used.
cplex.Model.lb = [-inf(1,wh+region_num) zeros(1,h_number)];
cplex.Model.ub = [inf(1,wh+region_num) inf(1,h_number)];
disp('done.');

fprintf('Initializing the constraint matrix...');
con_number = wh + 2*h_number; % The number of constraints.
% The number of non-zero elements in the DCT constraints between x and b.
NZ_number1 = wh2;
% The number of non-zero elements in the constraints between h and x.
NZ_number2 = 2*h_number*3;
% The maximal number of non-zero elements in the constraint matrix A.
NZ_number = NZ_number1 + NZ_number2;
% A vector of y-coordinates of all non-zero elements.
A_index_i = zeros(NZ_number,1);
% A vector of x-coordinates of all non-zero elements.
A_index_j = zeros(NZ_number,1);
% A vector of all non-zero elements.
A_value_ij = zeros(NZ_number,1);
cplex.Model.lhs = [img(:); -inf(2*h_number,1)];
cplex.Model.rhs = [img(:); zeros(2*h_number,1)];
% Initialize the number of constraints added.
NZ_index = 0;
con_ind = 0;
disp('done.');

% Add the constraints between x- and b-variables.
% x(i)=b(current region)+img(i)
% They are actually represented here as:
% img(i) <= x(i)-b(current region) <= img(i)
fprintf('Setting the constraints about differentially encoded DC coefficients...');
ri_matrix = blockIndexMatrix(h, w, regionSize_H, regionSize_W);
for j=1:w % Column by column
    for i=1:h % Row by row
        % Add the non-zero element corresponding to x(i).
        NZ_index = NZ_index + 1;
        con_ind = con_ind + 1;
        A_index_i(NZ_index) = con_ind;
        A_index_j(NZ_index) = con_ind; % = sub2ind([h w],i,j)
        A_value_ij(NZ_index) = 1;
        % Add the non-zero element corresponding to b(i).
        region_index = ri_matrix(i,j);
        NZ_index = NZ_index + 1;
        A_index_i(NZ_index) = con_ind;
        A_index_j(NZ_index) = wh + region_index;
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
h_ind = wh+region_num; % Initialize the index of the current h-variable
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
cplex.Model.ctype = repmat('C',[1 wh+region_num+h_number]);

Ts.seconds_model = toc(t1);

% Solve the optimization problem
disp('Solving the optimization problem...');
separator = repmat('=',1,80);
disp(separator);

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
objvals = [cplex_solution.objval sum(cplex_solution.x(wh+region_num+1:end)) im2smoothness(imgs.x)];
fprintf('Objective value of the returned solution: %g (Method: %s)\n', objvals(1), cplex_method);
fprintf('Objective value calculated from the retuned h-variables: %g\n', objvals(2));
fprintf('Objective value calculated from the retuned x-variables: %g\n\n', objvals(3));
assert(abs(objvals(1)-objvals(2))<eps && abs(objvals(2)-objvals(3))<eps, 'Objective values do not macth! CPLEX must have been running in an unstable mode! Results are not reliable!');

Ts.seconds_all = toc(t0);

disp(separator);

% Delete the log file generated by CPLEX as this is not needed.
% if exist('clone1.log','file')
%     delete('clone1.log');
% end
end % End of the main cplex_signbit() function

function bi_matrix = blockIndexMatrix(h,w,H,W,mark_mask)
% This internal function returns a matrix labeling the block index for each
% pixel in the blocks to be marked. If mark_mask is not specified, all
% blocks are to be numbered, if specified, only the blocks to be marked are
% counted, the indeces of the unmarked blocks are set to be 0

if ~exist('mark_mask','var')
    mark_mask = true(h/H,w/W);
end
bi_matrix = zeros(h,w);
block_index = 0;
block_index_marked = 0;
for j=1:W:w
    for i=1:H:h
        block_index = block_index + 1;
        if mark_mask(block_index)
            block_index_marked = block_index_marked + 1;
            bi_matrix(i:(i+H-1),j:(j+W-1)) = block_index_marked;
        end
    end
end
    
end % End of blockIndexMatrix() function