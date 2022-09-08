function plot_dctabs2dct(groud_truth, cplex_solution, mask)
% plot_dctabs2dct(groud_truth, cplex_solution, mask)
% This function plots the figure between the groud truth absolute values 
% and the recovered values (with signs) of 2-D DCT coefficients.
%
% Input arguments:
% groud_truth: ground truth value of DCT coefficients (this is needed
% rather than just the absolute values to know the groudth truth signs)
% cplex_solution: DCT values recovered by CPLEX
% mask: the mask used for recovering DC coefficients
%
% Shujun Li, Ruiyuan Lin @ University of Surrey, UK 2015
% Shujun Li: http://www.hooklee.com/

if nargin<2
    disp('Two input arguments are needed!');
    return;
end

if ~isequal(size(groud_truth),size(cplex_solution))
    disp('The two input arguments should be of the same size!');
    return;
end

h = size(groud_truth, 1);
w = size(groud_truth, 2);
mask = process_mask(mask, h, w);
if all(mask(:)) % All-true matrix
    disp('No any sign bits are unknown/recovered, so nothing can be shown!');
    return;
end

x = groud_truth(~mask(:));
x_abs = abs(x);
y = cplex_solution(~mask(:));
figure;
scatter(x_abs(x>0), y(x>0), '*b');
hold on;
scatter(x_abs(x<0), y(x<0), '+k');
scatter(x_abs(x==0), y(x==0), 'og');
x_abs_max = ceil(max(x_abs));
plot(0:x_abs_max, 0:x_abs_max, '--r');
plot(0:x_abs_max, 0:-1:-x_abs_max, '--r');
xlabel('Absolute (ground truth) value of DCT coefficients');
ylabel('Recovered values of DCT coeffients from CPLEX');
hold off;
axis equal;
axis([0 x_abs_max -x_abs_max x_abs_max]);

if exist('figureFullScreen.m', 'file')
    figureFullScreen(gcf);
end
