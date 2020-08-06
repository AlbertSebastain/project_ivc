function block = IAPBT(coeff)
%  Function Name : IDCT8x8.m
%  Input         : coeff (DCT Coefficients) 8*8*3
%  Output        : block (original image block) 8*8*3
[row,col,ded] = size(coeff);
iap = zeros(8,8);
for ii = 1:8
    for jj = 1:8
        i = ii-1;
        j = jj-1;
        if i == 0
            iapt = 1;
        else
            iapt = (2*8)/(8-j+sqrt(2)-1)*cos((j*(2*i+1)*pi)/(2*8));
        end
        iap(ii,jj) = iapt;
    end
end
for ii = 1 :ded
    block(:,:,ii) = iap*coeff(:,:,ii)*iap';
end
end