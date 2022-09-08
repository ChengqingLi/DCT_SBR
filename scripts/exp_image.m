% export restored images
stem = 'Lake';
U = 8;
DCTabs_min = 5;
DCpred_mode = 0;
[relaxX, relaxZ] = deal(false, false);
field = {'x1'}; % 'x', 'x00', 'x01'};
method = {'LP', 'naiveLP', 'negative', 'postive'};
for U = [1, 2, 3, 4, 5, 7]
    for DCpred_mode = 0:3
        file = dir(sprintf('../mat/lp/threshold/%s.pgm_*d%d_U%d_T%d_relaxX%d_relaxZ%d.mat',...
            stem, DCpred_mode, U, DCTabs_min, relaxX, relaxZ));
        data = load(fullfile(file.folder, file.name));
        for i = 1:length(field)
            suffix = field{i};
            name = sprintf('../images/output/%s_U%d_P%d_T%d_relaxX%d_relaxZ%d_%s.png', ...
                stem, U, DCpred_mode, DCTabs_min, relaxX, relaxZ, suffix);
            imwrite(data.imgs.(suffix), name);
            fprintf('[%s, U%d, P%d, %3s]: %.4f, %.2f\n', ...
                stem, U, DCpred_mode, method{i}, ...
                data.SSIMs.(suffix), data.PSNRs.(suffix));
        end
        fprintf('\n');
    end
end

% for DCpred_mode = 0:3
%     files = dir(sprintf('../mat/lp/threshold/*d%d_U%d_T%d_relaxX%d_relaxZ%d.mat',...
%         DCpred_mode, U, DCTabs_min, relaxX, relaxZ));
%     for i = 1:length(files)
%         if files(i).bytes < 4096
%             delete(fullfile(files(i).folder, files(i).name));
%         end
%     end
% end