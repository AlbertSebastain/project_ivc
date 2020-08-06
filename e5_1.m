%lena_small = double(imread('lena_small.tif'));
%% load image
for i = 20:40
    name = ['foreman00',num2str(i),'.bmp'];
    image(:,:,:,i-19) = double(imread(name));
    rgbimage(:,:,:,i-19) = image(:,:,:,i-19);
    image(:,:,:,i-19) = ictRGB2YCbCr(image(:,:,:,i-19));
    
end
%% initial parameters
[len,col,dem,fr] = size(image);
len_mv = [len/8*col/8];
%form       = double(imread('foreman0020.bmp'));
EOB = 4000; 
%scales = 1:0.6:1; % quantization scale factor, for E(4-1), we just evaluate scale factor of 1
%scales = [0.07, 0.2, 0.4, 0.8, 1.0, 1.5, 2, 3, 4, 4.5];
scales = 1;
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
    PSNR = zeros(1,21);
   % k_small  = IntraEncode(lena_small, qScale);
   %% encoding decoding first image
   fram_1 = image(:,:,:,1);
   [k_1,im_zigz1] = IntraEncode(fram_1,qScale,EOB);
   min_k(1) = min(k_1);
   len_k(1) = length(k_1);
   im_pmf = stats_marg(k_1,min(k_1):max(k_1));
   [Binarytree_first,HuffCode_first,BinCode_first,CodeLengths_first] = buildHuffman(im_pmf);
   bytestream{1} = enc_huffman_new(k_1-min(k_1)+1,BinCode_first,CodeLengths_first);
   k_rec_first = dec_huffman_new(bytestream{1},Binarytree_first,len_k(1));
   k_rec_first = k_rec_first+min_k(1)-1;
   im_rec(:,:,:,1) = IntraDecode(k_rec_first,size(fram_1),qScale,EOB);
   im_rec_rgb(:,:,:,1) = ictYCbCr2RGB(im_rec(:,:,:,1));
   min_mv = zeros(20);
   PSNR(1) = calcPSNR(rgbimage(:,:,:,1),im_rec_rgb(:,:,:,1));
   bits = numel(bytestream{1})*8;
   %bit_rate = bits/(len*col*dem/3)
   %% encoding and decoding follow images.
   for seq = 2:21

       fram = image(:,:,:,seq);
       ref_fram = im_rec(:,:,:,seq-1);
       fram_1dim = fram(:,:,1);
       ref_fram_1dim = ref_fram(:,:,1);
       mv = SSD(ref_fram_1dim,fram_1dim);
       res = fram-SSD_rec(ref_fram,mv); %the first res
       [res_k,im_zigz]    = IntraEncode(res, qScale,EOB);
       min_k(seq) = min(res_k);
       len_k(seq) = length(res_k);
        
  
       min_mv(seq-1) = min(mv(:));
       
       if seq == 2 % use the second image to construct huffman code.
           mv_pmf = stats_marg(mv,1:81);
           [Binarytree_mv,HuffCode_mv,BinCode_mv,CodeLengths_mv] = buildHuffman(mv_pmf);
           
           res_pmf = stats_marg(res_k,max_min:4000);
           [Binarytree_res,HuffCode_res,BinCode_res,CodeLengths_res] = buildHuffman(res_pmf);
       end
 
       %seq
       bytestream{seq} = enc_huffman_new(res_k-max_min+1,BinCode_res,CodeLengths_res); % huffman encode res.
       bytestream_mv{seq-1} = enc_huffman_new(mv(:)-min_mv(seq-1)+1,BinCode_mv,CodeLengths_mv); %huffman encode the vector
       k_rec = dec_huffman_new(bytestream{seq},Binarytree_res,len_k(seq));
       k_rec = k_rec+max_min-1;
       res_rec = IntraDecode(k_rec, size(fram_1),qScale,EOB);
       mv_rec = dec_huffman_new(bytestream_mv{seq-1},Binarytree_mv,len_mv);
       mv_rec = mv_rec+min_mv(seq-1)-1;
       mv_rec = reshape(mv_rec,[len/8,col/8]);
       im_rec_current = SSD_rec(ref_fram,mv_rec)+res_rec;
       im_rec(:,:,:,seq) = im_rec_current;
       im_rec_rgb(:,:,:,seq) = ictYCbCr2RGB(im_rec_current);
       PSNR(seq) = calcPSNR(rgbimage(:,:,:,seq),im_rec_rgb(:,:,:,seq));
       bits = bits+numel(bytestream{seq})*8+numel(bytestream_mv{seq-1})*8;
       %bit_rate = (numel(bytestream{seq})*8+numel(bytestream_mv{seq-1})*8)/(len*col)
       seq
   end
   PSNR_mean(scaleIdx) = mean(PSNR(:));
   bitPerPixel(scaleIdx) = bits/(numel(image)/3);
    fprintf('QP: %.1f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, bitPerPixel(scaleIdx), PSNR_mean(scaleIdx))
end
plot(bitPerPixel,PSNR_mean,'x-');
%% put all used sub-functions here.