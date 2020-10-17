function rec_image = SSD_rec(ref_image, motion_vectors)
%  Input         : ref_image(Reference Image, YCbCr image)
%                  motion_vectors
%
%  Output        : rec_image (Reconstructed current image, YCbCr image)
[len,col,dem] = size(ref_image);
rec_image = zeros(len,col,3);
ref_image = padarray(ref_image,[4,4],'both');
tx = len/8;
ty = col/8;
for i = 1:dem
    for j = 1:tx
        for k = 1:ty
            search = ref_image((j-1)*8+1:(j-1)*8+1+4+4+7,(k-1)*8+1:(k-1)*8+1+4+4+7,i);
            mv_ind = motion_vectors(j,k);
            mv = zeros(1,2);
            if mod(mv_ind,9) == 0
                mv(1) = mv_ind/9;
            else
                mv(1) = floor(mv_ind/9)+1;
            end
            mv(2) = mv_ind-(mv(1)-1)*9;
            search_block = search(mv(1):mv(1)+7,mv(2):mv(2)+7);
            block = search_block;
            rec_image((j-1)*8+1:j*8,(k-1)*8+1:k*8,i) = block;
        end
    end
end
end