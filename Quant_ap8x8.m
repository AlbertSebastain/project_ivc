function quant = Quant_ap8x8(dct_block, qScale)
%  Input         : dct_block (Original Coefficients, 8x8x3)
%                  qScale (Quantization Parameter, scalar)
%
%  Output        : quant (Quantized Coefficients, 8x8x3)
L = ones(8,8);
c = ones(8,8);
quant = zeros(8,8,3);
L = qScale.*L;
c = qScale.*c;
quant(:,:,1) = round(dct_block(:,:,1)./L);
quant(:,:,2) = round(dct_block(:,:,2)./c);
quant(:,:,3) = round(dct_block(:,:,3)./c);
end