function [blockwiseOrdered] = toBlockwiseOrder(twoDMatrix,H,W)
% [blockwiseOrdered] = toBlockwiseOrder(twoDMatrix,H,W)
% This function returns a column vector blockwiseOrdered containing
%     elements in twoDMatrix inserted in a blockwise manner (each block of
%     size [H W]): elements inside each block are inserted in a columnwise 
%     manner, blocks are scanned in a "columnwise" manner
%
% Ruiyuan Lin @ University of Surrey, UK 2015
[h, w] = size(twoDMatrix);
WH = W*H;
wh = w*h;
blockwiseOrdered = zeros(wh,1);
block_index = 0;
for j=1:W:w
   for i=1:H:h 
       block_index = block_index + 1;
       blockwiseOrdered((block_index-1)*WH + (1:WH)) = reshape(twoDMatrix(i:(i+H-1),j:(j+W-1)),[WH 1]);
   end
end

end