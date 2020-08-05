function rgb = ictYCbCr2RGB(yuv)
% Input         : yuv (Original YCbCr image)
% Output        : rgb (RGB Image after transformation)
% YOUR CODE HERE
Y = yuv(:,:,1);
Cb = yuv(:,:,2);
Cr = yuv(:,:,3);
r = Y+1.402.*Cr;
g = Y-0.344.*Cb-0.714.*Cr;
b = Y+1.772.*Cb;
rgb(:,:,1) = r;
rgb(:,:,2) = g;
rgb(:,:,3) = b;
end