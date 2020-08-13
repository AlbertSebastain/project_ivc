function quant = Quant8x8(dct_block, qScale,lum)
%  Input         : dct_block (Original Coefficients, 8x8x3)
%                  qScale (Quantization Parameter, scalar)
%                  lum "full" 3 dim quantization "luminance" luminance
%                  quantization "chorma" chorma quantization
%  Output        : quant (Quantized Coefficients, 8x8x3)

if nargin == 2
    lum = "full";
end
L = [16,11,10,16,24,40,51,61;12,12,14,19,26,58,60,55;14,13,16,24,40,57,69,56;14,17,22,29,...
    51,87,80,62;18,55,37,56,68,109,103,77;24,35,55,64,81,104,113,92;49,64,78,87,103,121,120,101;...
    72,92,95,98,112,100,103,99];
c = [17,18,24,47,99,99,99,99;18,21,26,66,99,99,99,99;24,13,56,99,99,99,99,99;47,66,99,99,99,99,99,99;99,99,...
    99,99,99,99,99,99;99,99,99,99,99,99,99,99;99,99,99,99,99,99,99,99;99,99,99,99,99,99,99,99];
[~,~,dim] = size(dct_block);
quant = zeros(8,8,dim);
L = qScale.*L;
c = qScale.*c;
% [~,~,dim] = size(dct_block);
if strcmp(lum,"full")
    quant_table = [L;c;c];
elseif strcmp(lum,"luminance")
    L = [1,1,1,1,2,4,5,6;1,1,1,1,2,5,6,5;1,1,1,2,4,5,6,5;1,1,2,2,...
    5,8,8,6;1,5,3,5,6,10,10,7;2,3,5,6,8,10,11,9;4,6,7,8,10,12,12,10;...
    7,9,9,9,11,10,10,9];
    L = L*qScale;
    quant_table = L;
else
    c = [1,1,2,4,9,9,9,9;1,2,2,6,9,9,9,9;2,1,5,9,9,9,9,9;4,6,9,9,9,9,9,9;9,9,...
    9,9,9,9,9,9;9,9,9,9,9,9,9,9;9,9,9,9,9,9,9,9;9,9,9,9,9,9,9,9];
    c = c*qScale;
    quant_table = c;
end
for d = 1:dim
    quant(:,:,d) = round(dct_block(:,:,d)./quant_table(d));
    %quant(:,:,2) = round(dct_block(:,:,2)./c);
    %quant(:,:,3) = round(dct_block(:,:,3)./c);
end