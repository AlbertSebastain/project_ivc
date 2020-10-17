function MSE = calcMSE(Image, recImage)
        [len,col,de] = size(Image);
    diff = Image - recImage;
    diff = diff.^2;
    sumd = sum(diff(:));
    MSE = sumd/(len*col*de);
end