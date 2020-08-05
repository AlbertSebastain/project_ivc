function coeff = DCT8x8(block)
%  Input         : block    (Original Image block, 8x8x3)
%
%  Output        : coeff    (DCT coefficients after transformation, 8x8x3)
dct8 = zeros(8,8);
[row,col,dem] = size(block);
coeff = zeros(row,col,dem);
for ii = 1:8
    for jj = 1:8
        k = ii-1;
        n = jj-1;
        if k == 0
            dct_t = 1/sqrt(8);
        else
            dct_t = sqrt(2/8)*cos(pi*(2*n+1)*k/(2*8));
        end
        dct8(ii,jj) = dct_t;
    end
end

for ii = 1:dem
    coeff(:,:,ii) = dct8*block(:,:,ii)*dct8.';
end
end