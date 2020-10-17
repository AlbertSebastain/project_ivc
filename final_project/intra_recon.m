function [rec] = intra_recon(img,im_size,qScale,EOB)
    len = im_size(1);
    col = im_size(2);
    dim = im_size(3);
    image_rund = ZeroRunDec_EoB(img,EOB);
    image_run_y = image_rund(1:len*col);
    image_run_cb = image_rund(len*col+1:len*col+len/2*col/2);
    image_run_cr = image_rund(len*col+len/2*col/2+1:end);
    image_run_y = reshape(image_run_y, [len/8*64,col/8]);
    image_run_cb = reshape(image_run_cb,[len/8/2*64,col/8/2]);
    image_run_cr = reshape(image_run_cr,[len/8*64/2,col/8/2]);
    image_dezigz_y = blockproc(image_run_y, [64,1],@(block_struct) DeZigZag8x8(block_struct.data));
    image_dezigz_cb = blockproc(image_run_cb, [64,1],@(block_struct) DeZigZag8x8(block_struct.data));
    image_dezigz_cr = blockproc(image_run_cr, [64,1],@(block_struct) DeZigZag8x8(block_struct.data));
    image_dequant_y = blockproc(image_dezigz_y,[8,8],@(block_struct) DeQuant8x8(block_struct.data,qScale,1,"luminance"));
    image_dequant_cb = blockproc(image_dezigz_cb,[8,8],@(block_struct) DeQuant8x8(block_struct.data,qScale,1,"chrom"));
    image_dequant_cr = blockproc(image_dezigz_cr,[8,8],@(block_struct) DeQuant8x8(block_struct.data,qScale,1,"chroma"));
    img_Rec_y = blockproc(image_dequant_y,[8,8],@(block_struct) IDCT8x8(block_struct.data));
    img_Rec_cb = blockproc(image_dequant_cb,[8,8],@(block_struct) IDCT8x8(block_struct.data));
    img_Rec_cr = blockproc(image_dequant_cr,[8,8],@(block_struct) IDCT8x8(block_struct.data));
    %img_Rec_y = image_dezigz_y;
    %img_Rec_cb = image_dezigz_cb;
    %img_Rec_cr = image_dezigz_cr;
    rec = zeros(im_size);
    para_y = [7/8,-1/2,5/8];
    para_cb = [3/8,-1/4,7/8];
    para_cr = para_cb;
    reco_y = dec_pre(img_Rec_y,para_y);
    reco_cb = dec_pre(img_Rec_cb,para_cb);
    reco_cr = dec_pre(img_Rec_cr,para_cr);
    reco_cb_u = resampling(reco_cb,[2,1]);
    reco_cr_u = resampling(reco_cr,[2,1]);
    rec(:,:,1) = reco_y;
    rec(:,:,2) = reco_cb_u;
    rec(:,:,3) = reco_cr_u;
end
function rec_im = dec_pre(res_im,para)
    [len,col] = size(res_im);
    rec_im = zeros(len,col);
    rec_im(1,:) = res_im(1,:);
    rec_im(:,1) = res_im(:,1);
    a1 = para(1);
    a2 = para(2);
    a3 = para(3);
    for i = 2:len
        for j = 2:col
            image_n = a1*rec_im(i,j-1)+a2*rec_im(i-1,j-1)+a3*rec_im(i-1,j);
            rec_im(i,j) = res_im(i,j)+image_n;
        end
    end
end
    