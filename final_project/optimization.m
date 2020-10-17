%lena_small = double(imread('lena_small.tif'));
clear
clc
%% load image
load('inter_intra_parameters.mat');
for i = 20:40
    name = ['foreman00',num2str(i),'.bmp'];
    image(:,:,:,i-19) = double(imread(name));
    rgbimage(:,:,:,i-19) = image(:,:,:,i-19);
    image(:,:,:,i-19) = ictRGB2YCbCr(image(:,:,:,i-19));
    
end
%% initial parameters
intra_mode = 1;
EOB = 4000; 
scales = [0.07, 0.2, 0.4, 0.8, 1.0, 1.5, 2, 3, 4, 4.5];
[len,col,dem,fr] = size(image);
len_mv = [len/8*col/8];
%form       = double(imread('foreman0020.bmp'));

 % quantization scale factor, for E(4-1), we just evaluate scale factor of 1

bitPerPixel = zeros(numel(scales),1);
PSNR_mean = zeros(numel(scales),1);
max_min = -1000;
%% initial coefficients
for scaleIdx = 1 : numel(scales)
   min_k = zeros(1,21);
   len_k = zeros(1,21);
    qScale   = scales(scaleIdx);
    bytestream = cell(21);
    bytestream_mv = cell(20);
    k = zeros(1,20);
    PSNR = zeros(1,21);
   % k_small  = IntraEncode(lena_small, qScale);
   %% encoding decoding first image
   fram_1 = image(:,:,:,1);
   if intra_mode == 1
        img_dst = intra_predict(fram_1,EOB,qScale);
   else
        [img_dst,im_zigz1] = IntraEncode(fram_1,qScale,EOB,intra_mode);
   end
   min_k(1) = min(img_dst);
   len_k(1) = length(img_dst);
   im_pmf = stats_marg(img_dst,min(img_dst):max(img_dst));
   [Binarytree_first,HuffCode_first,BinCode_first,CodeLengths_first] = buildHuffman(im_pmf);
   bytestream{1} = enc_huffman_new(img_dst-min(img_dst)+1,BinCode_first,CodeLengths_first);
   k_rec_first = dec_huffman_new(bytestream{1},Binarytree_first,len_k(1));
   k_rec_first = k_rec_first+min_k(1)-1;
   if intra_mode == 1
        im_rec(:,:,:,1) = intra_recon(k_rec_first,size(fram_1,[1,2,3]),qScale,EOB);
   else
        im_rec(:,:,:,1) = IntraDecode(k_rec_first,size(fram_1),qScale,EOB,intra_mode);
   end
   im_rec(:,:,:,1) = deblock(im_rec(:,:,:,1),1);
   im_rec_rgb(:,:,:,1) = ictYCbCr2RGB(im_rec(:,:,:,1));
   min_mv = zeros(20);
   PSNR(1) = calcPSNR(rgbimage(:,:,:,1),im_rec_rgb(:,:,:,1));
   bits = numel(bytestream{1})*8;
   bits_ebc = bits;
   %bit_rate = bits/(len*col*dem/3)
   %% encoding and decoding follow images.
   for seq = 2:21

       fram = image(:,:,:,seq);
       ref_fram = im_rec(:,:,:,seq-1);
       fram_1dim = fram(:,:,1);
       ref_fram_1dim = ref_fram(:,:,1);
       mv = SSD(ref_fram_1dim,fram_1dim);
       res = fram-SSD_rec(ref_fram,mv); %the first res
       [res_k,im_zigz]    = IntraEncode(res, qScale,EOB,intra_mode);
       min_k(seq) = min(res_k);
       len_k(seq) = length(res_k);
        
  
       min_mv(seq-1) = min(mv(:));
       
       if seq == 2 % use the second image to construct huffman code.
           %mv_pmf = stats_marg(mv,1:81);
           %[Binarytree_mv,HuffCode_mv,BinCode_mv,CodeLengths_mv] = buildHuffman(mv_pmf);
           res_pmf = stats_marg(res_k,max_min:4000);
           [Binarytree_res,HuffCode_res,BinCode_res,CodeLengths_res] = buildHuffman(res_pmf);
       end
 
       %seq
       bytestream{seq} = enc_huffman_new(res_k-max_min+1,BinCode_res,CodeLengths_res); % huffman encode res.
       %bytestream_mv{seq-1} = enc_huffman_new(mv(:)-min_mv(seq-1)+1,BinCode_mv,CodeLengths_mv); %huffman encode the vectort
       [bytestream_ebc{seq-1},bitb,k(seq-1)] = exponential_golomb_code(mv(:)-min_mv(seq-1)+1);
       k_rec = dec_huffman_new(bytestream{seq},Binarytree_res,len_k(seq));
       k_rec = k_rec+max_min-1;
       res_rec = IntraDecode(k_rec, size(fram_1),qScale,EOB,intra_mode);
       %mv_rec = dec_huffman_new(bytestream_mv{seq-1},Binarytree_mv,len_mv);
       mv_rec = exponential_golomb_decode(bytestream_ebc{seq-1},k(seq-1));
       mv_rec = mv_rec+min_mv(seq-1)-1;
       mv_rec = reshape(mv_rec,[len/8,col/8]);
       im_rec_current = SSD_rec(ref_fram,mv_rec)+res_rec;
       im_rec(:,:,:,seq) = im_rec_current;
       %im_rec_current = deblock(im_rec(:,:,:,seq),1);
       im_rec(:,:,:,seq) = deblock(im_rec(:,:,:,seq),1);
       im_rec_rgb(:,:,:,seq) = ictYCbCr2RGB(im_rec_current);
       PSNR(seq) = calcPSNR(rgbimage(:,:,:,seq),im_rec_rgb(:,:,:,seq));
       %bits = bits+numel(bytestream{seq})*8+numel(bytestream_mv{seq-1})*8;
       bits_ebc = bits_ebc+numel(bytestream{seq})*8+bitb;
       %bit_rate = (numel(bytestream{seq})*8+numel(bytestream_mv{seq-1})*8)/(len*col)
       seq
   end
   PSNR_mean(scaleIdx) = mean(PSNR(:));
   %bitPerPixel(scaleIdx) = bits/(numel(image)/3);
   bitPerPixel(scaleIdx) = bits_ebc/(numel(image)/3);
    fprintf('QP: %.1f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, bitPerPixel(scaleIdx), PSNR_mean(scaleIdx))
end
interval = bitPerPixel <4;
bitPerPixel_restrict = bitPerPixel(interval);
PSNR_mean_restrict = PSNR_mean(interval);
base_inter = bitrate_inter < 4;
bitrate_inter_restrict = bitrate_inter(base_inter);
PSNR_inter_restrict = PSNR_inter(base_inter);
interval_intra = bitrate_image < 4;
bitrate_intra_restrict = bitrate_image(interval_intra);
PSNR_intra_restrict = PSNR_image(interval_intra);
plot(bitPerPixel_restrict,PSNR_mean_restrict,'x-');
hold on
plot(bitrate_inter_restrict,PSNR_inter_restrict,'x-r');
hold on
plot(bitrate_intra_restrict,PSNR_intra_restrict,'x-g');
legend('optimization','baseline inter','baseline intra');
xlabel('bitrate');
ylabel('PSNR');
title('D R plot with restricted 0-4 bits');
hold off
figure
plot(bitPerPixel,PSNR_mean,'x-');
hold on
plot(bitrate_inter,PSNR_inter,'x-r');
plot(bitrate_image,PSNR_image,'x-g');
legend('optimization','baseline inter','baseline intra');
xlabel('bitrate');
ylabel('PSNR');
title('D R plot');
%% put all used sub-functions here.
