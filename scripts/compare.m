% MIP region size
DCpred_mode = 0;
for sz = [32]
    fprintf('%d*%d:\n', sz, sz);
    for U = [3, 5, 7]
        files = dir(sprintf('../mat/mip/region_size/*sz%d*%d*U%d*.mat', sz, sz, U));
        fileNum = length(files);
        [t, psnr, ssim] = deal(zeros(fileNum, 1));
        for k = 1:fileNum
            file = files(k);
            load(fullfile(file.folder, file.name));
            t(k) = Ts.seconds_1stPass + Ts.seconds_2ndPassMIP;
            psnr(k) = PSNRs.x1;
            ssim(k) = SSIMs.x1;
        end
        if fileNum ~= 0
            fprintf('[U%1d, P%1d, sz%d, %2d]: %8.3f;; ', U, DCpred_mode, sz, fileNum, mean(t));
            fprintf('%.6f, %.6f, %.4f, %.4f\n', mean(ssim), median(ssim), mean(psnr), median(psnr));
        end
    end
end

% LP
U = 50;
DCTabs_min = 5;
DCpred_mode = 0;
field = {'x00', 'x01', 'x', 'x1'};
method = {'negative', 'postive', 'naiveLP', 'LP'};
for U = [3, 5, 7]
    for i = 1:length(field)
        files = dir(sprintf('../mat/lp/threshold/*d%d_U%d_T%d_relaxX*_relaxZ*.mat',...
            DCpred_mode, U, DCTabs_min));
        fileNum = length(files);
        [psnr, ssim] = deal(zeros(fileNum, 1));
        for k = 1:fileNum
            file = files(k);
            load(fullfile(file.folder, file.name));
            [psnr(k), ssim(k)] = deal(PSNRs.(field{i}), SSIMs.(field{i}));
        end
        if fileNum ~= 0
            fprintf('[U%1d, P%1d, %2d, %2d, %3s]: ', ...
                U, DCpred_mode, DCTabs_min, fileNum, field{i});
            fprintf('%.6f, %.6f, %.4f, %.4f\n', mean(ssim), median(ssim), mean(psnr), median(psnr));
        end
    end
end
