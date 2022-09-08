DCTabs_min = 5;
for DCpred_mode = 1:2
    for U = [1, 2, 4, 8]
        for QF = [75, 85, 95]
            files = dir(sprintf('../mat/lp/Q%d/*d%d_U%d_T%d_relaxX1_relaxZ0.mat',...
                QF, DCpred_mode, U, DCTabs_min));
            fileNum = length(files);
            [psnr, ssim, t] = deal(zeros(fileNum, 1));
            for k = 1:fileNum
                file = files(k);
                load(fullfile(file.folder, file.name));
                t(k) = Ts.seconds_cplex;
                [psnr(k), ssim(k)] = deal(PSNRs.x1, SSIMs.x1);
            end
            if fileNum ~= 0
                fprintf('[%2d, P%1d, %2d, %2d, %2d]: ', ...
                    U, DCpred_mode, DCTabs_min, fileNum, QF);
                fprintf('%7.3f, %7.3f;; ', mean(t), max(t));
                fprintf('%.6f, %.6f, %7.4f, %7.4f\n', ...
                    mean(ssim), median(ssim), mean(psnr), median(psnr));
            end
        end
    end
end
