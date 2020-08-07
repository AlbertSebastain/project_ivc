function quant = Quant4x4(dct_block, qScale)
%  Input         : dct_block (Original Coefficients, 8x8x3)
%                  qScale (Quantization Parameter, scalar)
%
%  Output        : quant (Quantized Coefficients, 8x8x3)
L = [16,20,16,20;20,25,20,25;16,20,16,20;20,25,20,25];
c = [17,18,24,47;18,21,26,66;24,13,56,99;47,66,99,99];
quant = zeros(4,4,3);
L = qScale.*L;
c = qScale.*c;
quant(:,:,1) = round(dct_block(:,:,1)./L);
quant(:,:,2) = round(dct_block(:,:,2)./c);
quant(:,:,3) = round(dct_block(:,:,3)./c);
end