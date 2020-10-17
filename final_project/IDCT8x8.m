function block = IDCT8x8(coeff)
%  Function Name : IDCT8x8.m
%  Input         : coeff (DCT Coefficients) 8*8*3
%  Output        : block (original image block) 8*8*3
[row,col,ded] = size(coeff);
dct8 = zeros(8,8);
for i = 1:8
    for j = 1:8
        k = i-1;
        n = j-1;
        if k == 0
            dctte = 1/sqrt(8);
        else
            dctte = sqrt(2/8)*cos((2*n+1)*pi*k/(2*8));
        end
        dct8(i,j) = dctte;
    end
end
for i = 1 :ded
    block(:,:,i) = dct8'*coeff(:,:,i)*dct8;
end
end