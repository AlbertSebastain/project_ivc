function resam = resampling(mat,fac)
    mat1 = resample(mat,fac(1),fac(2),3);
    mat1 = resample(mat1',fac(1),fac(2),3);
    resam = mat1';
end