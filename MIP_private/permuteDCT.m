function [permuted_dct2] = permuteDCT(img_dct2, mask, universalPermute, H, W)
% [permuted_dct2] = permuteDCT(img_dct2, mask, universalPermute, H, W)
% This function permutes the input DCT coefficients blockwisely/universally
% 
% Input arguments:
%     img_dct2: input DCT coefficients
%     mask: mask of coefficients to be permuted, 0=to be permuted, 1=to remain
%     in original position
%     universalPermute: true: the permutation is done across the whole
%     image, false: the permutation is done within each block
%     H, W: size of each block
% Output arguments:
%     permuted_dct2: the permuted DCT coefficients
%
% Ruiyuan Lin @ University of Surrey, UK 2015
[h, w] = size(img_dct2);

permuted_dct2 = img_dct2;
if universalPermute
    permuteIndex = find(~mask);
    permuteList = zeros(numel(permuteIndex),1);
    permuteList(randperm(numel(permuteIndex))) = permuteIndex;
    permuted_dct2(permuteList) = img_dct2(permuteIndex);
else
    for j=1:W:w
        for i=1:H:h %blocks counted cloumnwise
            maskInBlock = mask(i:(i+H-1),j:(j+W-1));
            permuteIndexInBlock = find(~maskInBlock);
            if ~isempty(permuteIndexInBlock)
                permuteListInBlock = zeros(numel(permuteIndexInBlock),1);
                permuteListInBlock(randperm(numel(permuteIndexInBlock))) = permuteIndexInBlock; 
                img_dct2InBlock = img_dct2(i:(i+H-1),j:(j+W-1));
                permuted_dct2InBlock = img_dct2InBlock;
                permuted_dct2InBlock(permuteListInBlock) = img_dct2InBlock(permuteIndexInBlock);
                permuted_dct2(i:(i+H-1),j:(j+W-1)) = permuted_dct2InBlock;
            end
        end
    end
end

end