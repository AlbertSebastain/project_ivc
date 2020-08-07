function img = DCT4X4(im)
%     dctm4 = [1,1,1,1;2,1,-1,-2;1,-1,-1,1;1,-2,2,-1];
%     a = 1/2;
%     b = sqrt(2/5);
%     ef = [a^2,a*b/2,a^2,a*b/2;a*b/2,b^2/2,a*b/2,b^2;a^2,a*b,a^2,a*b;a*b/2,b^2/4,a*b/2,b^2/4];
img = zeros(4,4,3);
for ii = 1:4
    for jj = 1:4
        k = ii-1;
        n = jj-1;
        if k == 0
            dct_t = 1/sqrt(4);
        else
            dct_t = sqrt(2/4)*cos(pi*(2*n+1)*k/(2*4));
        end
        dct4(ii,jj) = dct_t;
    end
end
    for i = 1:3
        im_t = im(:,:,i);
    	img(:,:,i) = dct4*im_t*dct4';
    end
end
