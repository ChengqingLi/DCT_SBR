Us = [1, 4, 8];
% Us = [3, 5, 7, 9];

% DCTabs_min
U = 1;
for DCpred_mode = 0:3
    for DCTabs_min = 5
        files = dir(sprintf('../mat/lp/threshold/*d%d_U%d_T%d_relaxX*_relaxZ*.mat',...
        DCpred_mode, U, DCTabs_min));
        fileNum = length(files);
        [t, psnr, ssim, z_num] = deal(zeros(fileNum, 1));
        for k = 1:fileNum
            file = files(k);
            load(fullfile(file.folder, file.name));
            t(k) = Ts.seconds_cplex;
            z_num(k) = mask_neg_Z_num;
            [psnr(k), ssim(k)] = deal(PSNRs.x1, SSIMs.x1);
        end
        if fileNum ~= 0
            fprintf('[U%1d, P%1d, %2d, %2d]: %7.3f, %7.3f;; ', U, DCpred_mode, DCTabs_min, fileNum, mean(t), max(t));
            fprintf('%.6f, %.4f, %5.0f\n', mean(ssim), mean(psnr), mean(z_num));
        end
    end
end



% for i = 1:length(Us)
%     u = Us(i);
%     fprintf('%dU:\n', u);
%     for j = 1:length(SZs)
%         sz = SZs(j);
%         files = dir(sprintf('images/pgm/mip/*sz%d*%d*P1_D2*U%d*.mat', sz, sz, u));
%         [t1, t2, psnr, ssim] = deal(0);
%         fileNum = length(files);
%         for k = 1:fileNum
%             file = files(k);
%             load(fullfile(file.folder, file.name));
%             t1 = t1 + Ts.seconds_1stPass;
%             t2 = t2 + Ts.seconds_2ndPassMIP;
%             fields = {'all_max', 'all_maxname'};
%             psnr = psnr + struct2array(rmfield(PSNRs, fields));
%             ssim = ssim + struct2array(rmfield(SSIMs, fields));
%         end
%         if fileNum ~= 0
%             [t1, t2, psnr, ssim] = deal(t1/fileNum, t2/fileNum, psnr/fileNum, ssim/fileNum);
%             fprintf('%3d: %5.0f, %5.0f;;  ', sz, t1, t2);
%             fprintf('%.6f, %.6f, %.6f, %.6f;;  ', ssim(1), ssim(2), ssim(3), ssim(4));
%             fprintf('%.4f, %.4f, %.4f, %.4f\n', psnr(1), psnr(2), psnr(3), psnr(4));
%         end
%     end
% end

% % region size: first stage
% for i = 1:length(Us)
%     u = Us(i);
%     fprintf('%dU:\n', u);
%     for j = 1:length(SZs)
%         sz = SZs(j);
%         files = dir(sprintf('images/pgm/region_size/*sz%d*%d*U%d*.mat', sz, sz, u));
%         fileNum = length(files);
%         [t1, psnr, ssim] = deal(zeros(fileNum, 1));
%         for k = 1:fileNum
%             file = files(k);
%             load(fullfile(file.folder, file.name));
%             t1(k) = Ts.seconds_1stPass;
%             [psnr(k), ssim(k)] = deal(PSNRs.x, SSIMs.x);
%         end
%         if fileNum ~= 0
%             fprintf('%3d: %5.0f;;  ', sz, mean(t1));
%             fprintf('%.6f, %.6f;;  ', mean(ssim), median(ssim));
%             fprintf('%.4f, %.4f\n', mean(psnr), median(psnr));
%         end
%     end
% end

% for i = 1:length(Us)
%     u = Us(i);
%     fprintf('%dU:\n', u);
%     for j = 1:length(SZs)
%         sz = SZs(j);
%         files = dir(sprintf('images/pgm/mip/*sz%d*%d*U%d*.mat', sz, sz, u));
%         [t1, t2, psnr, ssim] = deal(0);
%         fileNum = length(files);
%         for k = 1:fileNum
%             file = files(k);
%             load(fullfile(file.folder, file.name));
%             t1 = t1 + Ts.seconds_1stPass;
%             t2 = t2 + Ts.seconds_2ndPassMIP;
%             fields = {'all_max', 'all_maxname'};
%             psnr = psnr + struct2array(rmfield(PSNRs, fields));
%             ssim = ssim + struct2array(rmfield(SSIMs, fields));
%         end
%         if fileNum ~= 0
%             [t1, t2, psnr, ssim] = deal(t1/fileNum, t2/fileNum, psnr/fileNum, ssim/fileNum);
%             fprintf('%3d: %5.0f, %5.0f;;  ', sz, t1, t2);
%             fprintf('%.6f, %.6f, %.6f, %.6f;;  ', ssim(1), ssim(2), ssim(3), ssim(4));
%             fprintf('%.4f, %.4f, %.4f, %.4f\n', psnr(1), psnr(2), psnr(3), psnr(4));
%         end
%     end
% end

% Print(5, 32)

% function  Print(U, sz)
%     for DCpred_mode = 1:3
%         fprintf('P%d:\n', DCpred_mode);
%         for DCdepend = 0:2
%             files = dir(sprintf('images/pgm/mip/*sz%d*P%d_D%d*U%d*.mat', sz, DCpred_mode, DCdepend, U));
%             [t1, t2, psnr, ssim] = deal(0);
%             fileNum = length(files);
%             for k = 1:fileNum
%                 file = files(k);
%                 load(fullfile(file.folder, file.name));
%                 t1 = t1 + Ts.seconds_1stPass;
%                 t2 = t2 + Ts.seconds_2ndPassMIP;
%                 fields = {'all_max', 'all_maxname'};
%                 psnr = psnr + struct2array(rmfield(PSNRs, fields));
%                 ssim = ssim + struct2array(rmfield(SSIMs, fields));
%             end
%             if fileNum ~= 0
%                 [t1, t2, psnr, ssim] = deal(t1/fileNum, t2/fileNum, psnr/fileNum, ssim/fileNum);
%                 fprintf('%1d: %5.0f, %5.0f;;  ', DCdepend, t1, t2);
%                 fprintf('%.6f, %.6f, %.6f, %.6f;;  ', ssim(1), ssim(2), ssim(3), ssim(4));
%                 fprintf('%.4f, %.4f, %.4f, %.4f\n', psnr(1), psnr(2), psnr(3), psnr(4));
%             end
%         end
%     end
% end

% DC coding
% for DCpred_mode = 1:3
%     fprintf('P%d:\n', DCpred_mode);
%     for DCdepend = 0:2
%         files = dir(sprintf('images/pgm/mip/*sz32*P%d_D%d*U5*.mat', DCpred_mode, DCdepend));
%         fileNum = length(files);
%         [t1, psnr, ssim] = deal(zeros(fileNum, 1));
%         for k = 1:fileNum
%             file = files(k);
%             load(fullfile(file.folder, file.name));
%             t1(k) = Ts.seconds_1stPass;
%             [psnr(k), ssim(k)] = deal(PSNRs.x1, SSIMs.x1);
%         end
%         if fileNum ~= 0
%             fprintf('%1d: %5.0f,;;  ', DCdepend, mean(t1));
%             fprintf('%.6f, %.6f;;  ', mean(ssim), median(ssim));
%             fprintf('%.4f, %.4f\n', mean(psnr), median(psnr));
%         end
%     end
% end
