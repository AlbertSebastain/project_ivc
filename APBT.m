function coeff = APBT(block)
%  Input         : block    (Original Image block, 8x8x3)
%
%  Output        : coeff    (DCT coefficients after transformation, 8x8x3)
ap8 = zeros(8,8);
[row,col,dem] = size(block);
coeff = zeros(row,col,dem);
for ii = 1:8
    for jj = 1:8
        i= ii-1;
        j = jj-1;
        if i == 0
            ap_t = 1/8;
        else
            ap_t = (8-i+sqrt(2)-1)/(8*8)*cos((i*(2*j+1)*pi)/(2*8));
        end
        ap8(ii,jj) = ap_t;
    end
end

for ii = 1:dem
    coeff(:,:,ii) = ap8*block(:,:,ii)*ap8.';
end
end