DCTabs_min = 5;
for DCpred_mode = 0:3
    for U = [1:8 10:5:30 40:10:50]
        files = dir(sprintf('../mat/lp/threshold/*d%d_U%d_T%d_relaxX0_relaxZ0.mat',...
            DCpred_mode, U, DCTabs_min));
        fileNum = length(files);
        [psnr, ssim, t] = deal(zeros(fileNum, 1));
        for k = 1:fileNum
            file = files(k);
            load(fullfile(file.folder, file.name));
            t(k) = Ts.seconds_cplex;
            [psnr(k), ssim(k)] = deal(PSNRs.x1, SSIMs.x1);
        end
        if fileNum ~= 0
            fprintf('[%2d, P%1d, %2d, %2d]: ', ...
                U, DCpred_mode, DCTabs_min, fileNum);
            fprintf('%7.3f, %7.3f;; ', mean(t), max(t));
            fprintf('%.6f, %.4f\n', mean(ssim), mean(psnr));
        end
    end
end
