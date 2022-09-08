% batch script to execute cplex_signbit.m
files = dir('images/pgm/*.pgm');
len = length(files);
[relaxX, relaxZ] = deal(false, false);
recovery_method = 0;

% global output_folder_name;
% output_folder_name = 'mat/lp/threshold';

for U = [1:4 5 6 7 8 10:5:30 40 50]
    for DCpred_mode = 0:3
        for DCTabs_min = 5
            for i = 1:len
                file = files(i);
                name = fullfile(file.folder, file.name);
                
                output = fullfile(file.folder, [file.name ...
                    '_signbitLP_DCpred' num2str(DCpred_mode) '_U' num2str(U) ...
                    '_T' num2str(DCTabs_min) '_relaxX' num2str(relaxX) ...
                    '_relaxZ' num2str(relaxZ) '.mat']);
                if exist(output, 'file')
                    fprintf('[U%d] skip %s, already existed\n', U, file.name);
                    continue;
                end

                mask = sprintf('masks/mask_%d.pgm', U);
                cplex_signbit(name, mask, DCpred_mode, DCTabs_min, ...
                    relaxX, relaxZ, recovery_method);
            end
        end
    end
end
% clear output_folder_name;
system('rm clone*.log');
