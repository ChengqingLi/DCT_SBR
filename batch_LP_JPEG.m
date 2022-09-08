% batch script to execute cplex_signbit.m
files = dir('images/pgm/*.pgm');
len = length(files);
DCpred_mode = 0;
DCTabs_min = 5;
[relaxX, relaxZ] = deal(false, false);
recovery_method = 0;

relax = [false, false; true, false; false, true; true, true];
QF = 95;
addpath('./jpeg_toolbox');
% global output_folder_name;
% output_folder_name = sprintf('mat/lp/Q%d', QF);

for U = [1:4 5 6 7 8 10:5:30 40 50]
    for relax_index = 2
        [relaxX, relaxZ] = deal(relax(relax_index, 1), relax(relax_index, 2));
        for DCpred_mode = 1:2
            for i = 1:len
                file = files(i);
                [~, stem] = fileparts(file.name);
                
                image = imread(fullfile(file.folder, file.name));
                jpeg_folder = sprintf('images/JPEG/Q%d/', QF);
                jpeg_name = sprintf('%s/%s.jpg', jpeg_folder, stem);
                imwrite(image, jpeg_name, 'Quality', QF);
      
                output = fullfile(jpeg_folder, [stem '.jpg' ...
                    '_signbitLP_DCpred' num2str(DCpred_mode) '_U' num2str(U) ...
                    '_T' num2str(DCTabs_min) '_relaxX' num2str(relaxX) ...
                    '_relaxZ' num2str(relaxZ) '.mat']);
                if exist(output, 'file')
                    fprintf('[U%d] skip %s, already existed\n', U, file.name);
                    continue;
                end

                mask = sprintf('masks/mask_%d.pgm', U);
                cplex_signbit(jpeg_name, mask, DCpred_mode, DCTabs_min, ...
                    relaxX, relaxZ, recovery_method);
            end
        end
    end
end
% clear output_folder_name;
system('rm clone*.log');
