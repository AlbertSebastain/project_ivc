function im = IDCT4x4(img)
    im = zeros(4,4,3);
%     idct4 = [1,1,1,1/2;1,1/2,-1,-1;1,-1/2,-1,1;1,-1,1,-1/2];
%     a = 1/2;
%     b = sqrt(2/5);
%     ef = [a^2,a*b,a^2,a*b;a*b,b^2,a*b,b^2;a^2,a*b,a^2,a*b;a*b,b^2,a*b,b^2];
for i = 1:4
    for j = 1:4
        k = i-1;
        n = j-1;
        if k == 0
            dctte = 1/sqrt(4);
        else
            dctte = sqrt(2/4)*cos((2*n+1)*pi*k/(2*4));
        end
        dct4(i,j) = dctte;
    end
end
    for i = 1:3
        img_t = img(:,:,i);
        im(:,:,i) = dct4'*img_t*dct4;
    end
end