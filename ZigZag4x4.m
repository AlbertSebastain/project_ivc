function zz = ZigZag4x4(quant)
%  Input         : quant (Quantized Coefficients, 8x8x3)
%
%  Output        : zz (zig-zag scaned Coefficients, 64x3)
zigz =     [1,2,6,7;
            3,5,8,13;
            4,9,12,14;
            10,11,15,16];
for k = 1:3
    quantt = quant(:,:,k);
    zz(zigz(:),k) = quantt(:);
end
end