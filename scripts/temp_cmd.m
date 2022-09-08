files = dir('../mat/mip/region_size/*.mat');
for i = 1:length(files)
    file = files(i);
    if mod(i, 5) == 0 && file.bytes < 1024 * 1024 * 10
        continue;
    end
    name = fullfile(file.folder, file.name);
    load(name);
    save(name, 'imgs', 'Ts', 'SB', 'PSNRs', 'SSIMs', 'BERs', ...
    'DCpred_mode', 'DCTabs_min', 'recovery_method', 'DC_shift_level', 'H', 'W', 'Xmax', 'timeLimit', ...
    'TsSeg', 'relaxXBound', 'relaxXBound2');
end