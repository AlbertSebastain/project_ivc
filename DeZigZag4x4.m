function coeffs = DeZigZag4x4(zz)
%  Function Name : DeZigZag8x8.m
%  Input         : zz    (Coefficients in zig-zag order)
%
%  Output        : coeffs(DCT coefficients in original order)
zigz =     [1,2,6,7;
            3,5,8,13;
            4,9,12,14;
            10,11,15,16];
for k = 1:3
    coeffs_t = zz(zigz(:),k);
    coeffs_t = reshape(coeffs_t,[4,4]);
    coeffs(:,:,k) = coeffs_t;
end
end