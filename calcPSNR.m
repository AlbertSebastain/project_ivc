function PSNR = calcPSNR(Image, recImage)
    MSE = calcMSE(Image,recImage);
    PSNR = 10*log10((2^8-1)^2/MSE);
end