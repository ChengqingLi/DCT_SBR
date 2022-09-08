U = 2;
DCTabs_min = 5;
for DCpred_mode = 0:3
    files = dir(sprintf('../mat/lp/threshold/*d%d_U%d_T%d_relaxX*_relaxZ*.mat',...
        DCpred_mode, U, DCTabs_min));
    fileNum = length(files);
    [psnr, ssim] = deal(zeros(fileNum, 1));
    for k = 1:fileNum
        file = files(k);
        load(fullfile(file.folder, file.name));
        [psnr(k), ssim(k)] = deal(PSNRs.x1, SSIMs.x1);
    end
    if fileNum ~= 0
        fprintf('[U%1d, P%1d, %2d, %2d]: ', U, DCpred_mode, DCTabs_min, fileNum);
        fprintf('%.6f, %.6f, %.4f, %.4f\n', mean(ssim), median(ssim), mean(psnr), median(psnr));
    end
end