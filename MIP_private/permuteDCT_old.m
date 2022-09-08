function [permuted_dct2, permuteMatrix, invPermuteMatrix] = permuteDCT(img_dct2, mask, universalPermute, H, W)
% [permuted_dct2, permuteMatrix, invPermuteMatrix] = permuteDCT(img_dct2, mask, universalPermute, H, W)
% This function permutes the input DCT coefficients blockwisely/universally
% 
% Input arguments:
%     img_dct2: input DCT coefficients
%     mask: mask of coefficients to be permuted, 0=to be permuted, 1=to remain
%     in original position
%     H, W: size of each block
% Output arguments:
%     permuted_dct2: the permuted DCT coefficients
%     permuteMatrix: the permutation matrix the permutation is based on 
%     (i.e. permuted_dct2 = permuteMatrix * img_dct2)
%     when universalPermute = false, each block has its own 
%     permuteMatrixInBlock, the output permuteMatrix of size [WH wh] is 
%     formed by concatenating all the permuteMatrixInBlock horizontally, 
%     blocks are ordered "columnwisely"
%     invPermuteMatrix: the inverse of permutation matrix of each block
%     concatenated in the same manner as permuteMatrix
%
% Ruiyuan Lin @ University of Surrey, UK 2015
[h, w] = size(img_dct2);
WH = W*H;
wh = w*h;
permuted_dct2 = img_dct2;
if universalPermute
    permuteIndex = find(mask==0);
    permuteMatrix = eye(wh);
    permuteList = zeros(numel(permuteIndex),1);
    permuteList(randperm(numel(permuteIndex))) = permuteIndex;
    permuteMatrix(permuteList,:) = permuteMatrix(permuteIndex,:);
    permuted_dct2(:) = permuteMatrix * img_dct2(:);
    invPermuteMatrix = permuteMatrix'; %permutation matrix is non-signular, PP'=I
else
    block_index = 0;
    permuteMatrix = zeros(WH,wh);
    invPermuteMatrix = zeros(WH,wh);
    for j=1:W:w
        for i=1:H:h %blocks counted cloumnwise
            block_index = block_index + 1;
            maskInBlock = mask(i:(i+H-1),j:(j+W-1));
            permuteIndexInBlock = find(maskInBlock==0);
            if ~isempty(permuteIndexInBlock)
                permuteListInBlock = zeros(numel(permuteIndexInBlock),1);
                permuteListInBlock(randperm(numel(permuteIndexInBlock))) = permuteIndexInBlock; 
                permuteMatrixInBlock = eye(WH);
                permuteMatrixInBlock(permuteListInBlock,:) = permuteMatrixInBlock(permuteIndexInBlock,:);
                permuteMatrix(:,(block_index-1)*WH+(1:WH)) = permuteMatrixInBlock;
                invPermuteMatrix(:,(block_index-1)*WH+(1:WH)) = permuteMatrixInBlock';
                img_dct2InBlock = img_dct2(i:(i+H-1),j:(j+W-1));
                permuted_dct2InBlock = img_dct2InBlock;
                permuted_dct2InBlock(:) = permuteMatrixInBlock * img_dct2InBlock(:);
                permuted_dct2(i:(i+H-1),j:(j+W-1)) = permuted_dct2InBlock;
            else
                permuteMatrix(:,(block_index-1)*WH+(1:WH)) = eye(WH);
                invPermuteMatrix(:,(block_index-1)*WH+(1:WH)) = eye(WH);
            end
        end
    end
end
end