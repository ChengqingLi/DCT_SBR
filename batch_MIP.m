% batch script to execute cplex_signbitHierarchyNoNaiveLP.m
files = dir('images/pgm/*.pgm');
len = length(files);
DCpred_mode = 2;
DCdepends = [0, 1, 2];
DCTabs_min = 0;
[relaxXBound, relaxZBound] = deal(true, false);
DCpassByRegion = 0;  %  1:region-based, 2:block-based
timeLimit = 600;
relatGap = [];
DCpass_method = 0;  % 1:MIP, 2:LP
recovery_method = 0;

% global output_folder_name;
% output_folder_name = 'mat/mip/region_size';

for sz = [16, 32, 64, 128]
    [regionSize_H, regionSize_W] = deal(sz);
    for U = [3, 5, 7]
        for DCdepends_index = 1:length(DCdepends)
            DCdepend = DCdepends(DCdepends_index);
            for i = 16:len 
                file = files(i);
                name = fullfile(file.folder, file.name);
                
                output = fullfile(file.folder, [file.name '_signbitHMIP_sz' ...
                        num2str(regionSize_H) 'x' num2str(regionSize_W) '_P' num2str(DCpred_mode) ...
                        '_D' num2str(DCdepend) '_U' num2str(U) '_T' num2str(DCTabs_min) ...
                        '_tL' num2str(timeLimit) '_gap' num2str(relatGap) ...
                        '_rX' num2str(relaxXBound) '_rZ' num2str(relaxZBound) ...
                        '_DCpassByRegion' num2str(DCpassByRegion) '.mat']);

                if exist(output, 'file')
                    fprintf('[%dx%d U%d] skip %s, already existed\n', ...
                        regionSize_H, regionSize_W, U, file.name);
                    continue;
                end
                mask = sprintf('masks/mask_%d.pgm', U);
                [imgs, Ts, imgs_dct2, PSNRs, SSIMs, BERs] = cplex_signbitHierarchyDCACMIP(...
                    name, mask, regionSize_H, regionSize_W, DCpred_mode, DCdepend, ...
                    relaxXBound, relaxZBound, DCpassByRegion, DCTabs_min, timeLimit, relatGap, ...
                    DCpass_method, recovery_method);
            end
        end
    end
end
% clear output_folder_name;
system('rm clone*.log');