function dct_block = DeQuant_ap8x8(quant_block, qScale)
%  Function Name : DeQuant8x8.m
%  Input         : quant_block  (Quantized Block, 8x8x3)
%                  qScale       (Quantization Parameter, scalar)
%
%  Output        : dct_block    (Dequantized DCT coefficients, 8x8x3)
L = ones(8,8);
c = ones(8,8);
L = qScale.*L;
c = qScale.*c;
dct_block(:,:,1) = quant_block(:,:,1).*L;
dct_block(:,:,2) = quant_block(:,:,2).*c;
dct_block(:,:,3) = quant_block(:,:,3).*c;
end