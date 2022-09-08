U = 2;
DCTabs_min = 5;
relaxs = [false, false; true, false; false, true; true, true];
for DCpred_mode = 1:2
    for relax_index = 1:length(relaxs)
        [relaxX, relaxZ] = deal(relaxs(relax_index, 1), relaxs(relax_index, 2));
        files = dir(sprintf('../mat/lp/relax/*d%d_U%d_T%d_relaxX%d_relaxZ%d.mat',...
            DCpred_mode, U, DCTabs_min, relaxX, relaxZ));
        fileNum = length(files);
        [psnr, ssim, t] = deal(zeros(fileNum, 1));
        for k = 1:fileNum
            file = files(k);
            load(fullfile(file.folder, file.name));
            t(k) = Ts.seconds_cplex;
            [psnr(k), ssim(k)] = deal(PSNRs.x1, SSIMs.x1);
        end
        if fileNum ~= 0
            fprintf('[U%1d, P%1d, %2d, %2d, %d, %d]: ', ...
                U, DCpred_mode, DCTabs_min, fileNum, relaxX, relaxZ);
            fprintf('%7.3f, %7.3f;; ', median(t), max(t));
            fprintf('%.6f, %.6f, %.4f, %.4f\n', mean(ssim), median(ssim), mean(psnr), median(psnr));
        end
    end
end
