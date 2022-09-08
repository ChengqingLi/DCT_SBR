% MIP region size
for sz = [8, 16, 32, 64]
    fprintf('%d*%d:\n', sz, sz);
    for U = [3, 5, 7, 9]
        files = dir(sprintf('../mat/mip/region_size/*sz%d*%d*U%d*.mat', sz, sz, U));
        [t1, t2, psnr, ssim] = deal(0);
        fileNum = length(files);
        for k = 1:fileNum
            file = files(k);
            load(fullfile(file.folder, file.name));
            t1 = t1 + Ts.seconds_1stPass;
            t2 = t2 + Ts.seconds_2ndPassMIP;
            fields = {'all_max', 'all_maxname'};
            psnr = psnr + struct2array(rmfield(PSNRs, fields));
            ssim = ssim + struct2array(rmfield(SSIMs, fields));
        end
        if fileNum ~= 0
            [t1, t2, psnr, ssim] = deal(t1/fileNum, t2/fileNum, psnr/fileNum, ssim/fileNum);
            fprintf('%dU: %5.0f, %5.0f;;  ', U, t1, t2);
            fprintf('%.3f, %.3f, %.3f, %.3f, %.3f;;  ', ssim(1), ssim(2), ssim(3), ssim(4), max(ssim));
            fprintf('%.3f, %.3f, %.3f, %.3f, %.3f\n', psnr(1), psnr(2), psnr(3), psnr(4), max(psnr));
        end
    end
end