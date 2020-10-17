function motion_vectors_indices = SSD(ref_image, image)
%  Input         : ref_image(Reference Image, size: height x width)
%                  image (Current Image, size: height x width)
%
%  Output        : motion_vectors_indices (Motion Vector Indices, size: (height/8) x (width/8) x 1 )ref
range = 4;
ref_image = padarray(ref_image,[range,range],0,'both');
%ref_image = padarray(ref_image,[7,7],0,'post');
[len,col] = size(image);
tx = len/8;
ty = col/8;
for i = 1:tx
    for j = 1:ty
        block = image((i-1)*8+1:i*8,(j-1)*8+1:j*8);
        search = ref_image((i-1)*8+1:(i-1)*8+1+range+range+7,(j-1)*8+1:(j-1)*8+1+range+range+7);
        mv = searchmv(search,block);
        mv_ind = (mv(1)-1+5)*9+mv(2)+5;
        motion_vectors_indices(i,j) = mv_ind;
    end
end
end