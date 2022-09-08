function zz = zigzagMatrix(H,W)
% zz = zigzagMatrix(H,W)
% This function returns a zigzag traversal matrix (0-based) of size [H W]
% 
% Ruiyuan Lin @ University of Surrey, UK 2015
    WH = W*H;
    v=zigzag(reshape((1:WH),[H W]));
    zz=zeros(H,W);
    zz(v)=0:(WH-1);
end