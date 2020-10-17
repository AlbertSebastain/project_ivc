function pmf = stats_marg(image,range) 
    [cc] = hist(image(:),range);
    pmf = cc/sum(cc);

end