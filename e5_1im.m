%lena_small = double(imread('lena_small.tif'));
Lena = double(imread('lena.tif'));
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
scales = [0.07, 0.2, 0.4, 0.8, 1.0, 1.5, 2, 3, 4, 4.5];
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
   bitrate_im = zeros(21,1);
   PSNR_im = zeros(21,1);
   bitrate_im(1) = bits/(len*col);
   PSNR_im(1) = PSNR(1);
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
       [bitrate_im(seq),PSNR_im(seq)] = im_code(qScale,fram,EOB);
       seq
   end
   PSNR_image(scaleIdx) = sum(PSNR_im)/21;
   bitrate_image(scaleIdx) = sum(bitrate_im)/21;
   PSNR_mean(scaleIdx) = mean(PSNR(:));
   bitPerPixel(scaleIdx) = bits/(numel(image)/3);
    fprintf('QP: %.1f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, bitPerPixel(scaleIdx), PSNR_mean(scaleIdx))
end
plot(bitPerPixel,PSNR_mean,'x-');
hold on
title('dr plot');
xlabel('bitPerPixel');
ylabel('PSNR');
plot(bitrate_image,PSNR_image,'x-');

%% put all used sub-functions here.
function [Bitrate,PSNR_IM] = im_code(qScale,Lena,EOB)
[k,im_zigz]    = IntraEncode(Lena, qScale,EOB);
    im_pmf = stats_marg(k,max(k)-min(k)+1);
    [Binarytree,HuffCode,BinCode,CodeLengths] = buildHuffman(im_pmf);
    bytestream = enc_huffman_new(k-min(k)+1,BinCode,CodeLengths);
    %% use pmf of k_small to build and train huffman table
    %your code here
    %% use trained table to encode k to get the bytestream
    % your code here

    k_rec = dec_huffman_new(bytestream,Binarytree,length(k));
    k_rec = k_rec+min(k)-1;
    Bitrate = (numel(bytestream)*8) / (numel(Lena)/3);
    %% image reconstruction
    I_rec = IntraDecode(k_rec, size(Lena),qScale,EOB);
    I_rec = ictYCbCr2RGB(I_rec);
    Lena = ictYCbCr2RGB(Lena);
    PSNR_IM = calcPSNR(I_rec,Lena);
    %PSNR(scaleIdx) = calcPSNR(Lena, I_rec);
    %fprintf('QP: %.1f bit-rate: %.2f bits/pixel PSNR: %.2fdB\n', qScale, bitPerPixel(scaleIdx), PSNR(scaleIdx))
end
    function dst = IntraDecode(image, img_size , qScale,EOB)
%  Function Name : IntraDecode.m
%  Input         : image (zero-run encoded image, 1xN)
%                  img_size (original image size)
%                  qScale(quantization scale)
%  Output        : dst   (decoded image)
    image_rund = ZeroRunDec_EoB(image,EOB);
    image_rund = reshape(image_rund, [img_size(1)/8*64,img_size(2)/8*img_size(3)]);
    image_dezigz = blockproc(image_rund, [64,3],@(block_struct) DeZigZag8x8(block_struct.data));
    image_dequant = blockproc(image_dezigz,[8,8],@(block_struct) DeQuant8x8(block_struct.data,qScale));
    image_idc = blockproc(image_dequant,[8,8],@(block_struct) IDCT8x8(block_struct.data));
    %dst = ictYCbCr2RGB(image_idc);
    dst = image_idc;
end

function [dst,im_zigz] = IntraEncode(image, qScale,EOB)
%  Function Name : IntraEncode.m
%  Input         : image (Original RGB Image)
%                  qScale(quantization scale)
%  Output        : dst   (sequences after zero-run encoding, 1xN)
%im_tran = ictRGB2YCbCr(image);
im_tran = image;
[row,col,dem] = size(im_tran);
im_dct = blockproc(im_tran,[8,8],@(block_struct) DCT8x8(block_struct.data));
im_quant = blockproc(im_dct,[8,8],@(block_struct) Quant8x8(block_struct.data,qScale));
im_zigz = blockproc(im_quant,[8,8],@(block_struct) ZigZag8x8(block_struct.data));
im_zigz = im_zigz(:);
dst = ZeroRunEnc_EoB(im_zigz,EOB);
end

%% and many more functions

function coeff = DCT8x8(block)
%  Input         : block    (Original Image block, 8x8x3)
%
%  Output        : coeff    (DCT coefficients after transformation, 8x8x3)
dct8 = zeros(8,8);
[row,col,dem] = size(block);
coeff = zeros(row,col,dem);
for ii = 1:8
    for jj = 1:8
        k = ii-1;
        n = jj-1;
        if k == 0
            dct_t = 1/sqrt(8);
        else
            dct_t = sqrt(2/8)*cos(pi*(2*n+1)*k/(2*8));
        end
        dct8(ii,jj) = dct_t;
    end
end

for ii = 1:dem
    coeff(:,:,ii) = dct8*block(:,:,ii)*dct8.';
end
end
function block = IDCT8x8(coeff)
%  Function Name : IDCT8x8.m
%  Input         : coeff (DCT Coefficients) 8*8*3
%  Output        : block (original image block) 8*8*3
[row,col,ded] = size(coeff);
dct8 = zeros(8,8);
for i = 1:8
    for j = 1:8
        k = i-1;
        n = j-1;
        if k == 0
            dctte = 1/sqrt(8);
        else
            dctte = sqrt(2/8)*cos((2*n+1)*pi*k/(2*8));
        end
        dct8(i,j) = dctte;
    end
end
for i = 1 :ded
    block(:,:,i) = dct8'*coeff(:,:,i)*dct8;
end
end
function zze = ZeroRunEnc_EoB(zz,EOB)
%EOB = 1000;
%  Input         : zz (Zig-zag scanned sequence, 1xN)
%                  EOB (End Of Block symbol, scalar)
%
%  Output        : zze (zero-run-level encoded sequence, 1xM)
len = length(zz);
zze = [];
for rec = 1:len/64
    if rec ~= len/64
        zz_t = zz((rec-1)*64+1:rec*64);
    else
        zz_t = zz((rec-1)*64+1:len);
    end

[lens] = length(zz_t);
k = 1;
zze_t = [];
flag = 0;
m = 0;
while k <= lens
    if zz_t(k) == 0
        count = -1;
        m = k;
        while (zz_t(m) == 0 )
            if m == 1
                flag1 = 1;
            end
            count = count+1;
            m = m+1;
            if m > lens
                flag = 1;
                break
            end
            
        end
        zze_t = [zze_t,0,count];
        k = m-1;
    else
        zze_t = [zze_t,zz_t(k)];
    end
    k = k+1;
end
len2 = length(zze_t);
if flag == 1
    zze_t = zze_t(1:len2-1);
    zze_t(end) = EOB;
end
zze = [zze,zze_t];
end
end
function dst = ZeroRunDec_EoB(src,EOB)
%EoB = 1000;
%  Function Name : ZeroRunDec1.m zero run level decoder
%  Input         : src (zero run encoded sequence 1xM with EoB signs)
%                  EoB (end of block sign)
%
%  Output        : dst (reconstructed zig-zag scanned sequence 1xN)
k = 0;
len_src = length(src);
m = 1;
dst = [];
while k < len_src
    k = k+1;
    if m == 65
        m = 1;
    end
    if src(k) == EOB
        dst(end+1:end+64-m+1) = 0;
        m = 1;
        continue
    end
    if src(k) == 0
        dst(end+1:end+src(k+1)+1) = 0;
        m = m+src(k+1)+1;
        k = k+1;
       
        continue
    end
    dst(end+1) = src(k);
    m = m+1;
end
size(dst);
size(src);
end
function pmf = stats_marg(image,range) 
    [cc] = hist(image(:),range);
    pmf = cc/sum(cc);

end
function yuv = ictRGB2YCbCr(rgb)
% Input         : rgb (Original RGB Image)
% Output        : yuv (YCbCr image after transformation)
% YOUR CODE HERE
yuv(:,:,1) = 0.299.*rgb(:,:,1)+0.587.*rgb(:,:,2)+0.114.*rgb(:,:,3);
yuv(:,:,2) = -0.169.*rgb(:,:,1)-0.331.*rgb(:,:,2)+0.5.*rgb(:,:,3);
yuv(:,:,3) = 0.5.*rgb(:,:,1)-0.419.*rgb(:,:,2)-0.081.*rgb(:,:,3);
end
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
function PSNR = calcPSNR(Image, recImage)
    MSE = calcMSE(Image,recImage);
    PSNR = 10*log10((2^8-1)^2/MSE);
end

function MSE = calcMSE(Image, recImage)
        [len,col,de] = size(Image);
    diff = Image - recImage;
    diff = diff.^2;
    sumd = sum(diff(:));
    MSE = sumd/(len*col*de);
end
function quant = Quant8x8(dct_block, qScale)
%  Input         : dct_block (Original Coefficients, 8x8x3)
%                  qScale (Quantization Parameter, scalar)
%
%  Output        : quant (Quantized Coefficients, 8x8x3)
L = [16,11,10,16,24,40,51,61;12,12,14,19,26,58,60,55;14,13,16,24,40,57,69,56;14,17,22,29,...
    51,87,80,62;18,55,37,56,68,109,103,77;24,35,55,64,81,104,113,92;49,64,78,87,103,121,120,101;...
    72,92,95,98,112,100,103,99];
c = [17,18,24,47,99,99,99,99;18,21,26,66,99,99,99,99;24,13,56,99,99,99,99,99;47,66,99,99,99,99,99,99;99,99,...
    99,99,99,99,99,99;99,99,99,99,99,99,99,99;99,99,99,99,99,99,99,99;99,99,99,99,99,99,99,99];
quant = zeros(8,8,3);
L = qScale.*L;
c = qScale.*c;
quant(:,:,1) = round(dct_block(:,:,1)./L);
quant(:,:,2) = round(dct_block(:,:,2)./c);
quant(:,:,3) = round(dct_block(:,:,3)./c);
end
function dct_block = DeQuant8x8(quant_block, qScale)
%  Function Name : DeQuant8x8.m
%  Input         : quant_block  (Quantized Block, 8x8x3)
%                  qScale       (Quantization Parameter, scalar)
%
%  Output        : dct_block    (Dequantized DCT coefficients, 8x8x3)
L = [16,11,10,16,24,40,51,61;12,12,14,19,26,58,60,55;14,13,16,24,40,57,69,56;14,17,22,29,...
    51,87,80,62;18,55,37,56,68,109,103,77;24,35,55,64,81,104,113,92;49,64,78,87,103,121,120,101;...
    72,92,95,98,112,100,103,99];
c = [17,18,24,47,99,99,99,99;18,21,26,66,99,99,99,99;24,13,56,99,99,99,99,99;47,66,99,99,99,99,99,99;99,99,...
    99,99,99,99,99,99;99,99,99,99,99,99,99,99;99,99,99,99,99,99,99,99;99,99,99,99,99,99,99,99];
L = qScale.*L;
c = qScale.*c;
dct_block(:,:,1) = quant_block(:,:,1).*L;
dct_block(:,:,2) = quant_block(:,:,2).*c;
dct_block(:,:,3) = quant_block(:,:,3).*c;
end
function zz = ZigZag8x8(quant)
%  Input         : quant (Quantized Coefficients, 8x8x3)
%
%  Output        : zz (zig-zag scaned Coefficients, 64x3)
zigz =     [1     2     6     7    15    16    28    29;
     3     5     8    14    17    27    30    43;
     4     9    13    18    26    31    42    44;
    10    12    19    25    32    41    45    54;
    11    20    24    33    40    46    53    55;
    21    23    34    39    47    52    56    61;
    22    35    38    48    51    57    60    62;
    36    37    49    50    58    59    63    64];
for k = 1:3
    quantt = quant(:,:,k);
    zz(zigz(:),k) = quantt(:);
end
end
function coeffs = DeZigZag8x8(zz)
%  Function Name : DeZigZag8x8.m
%  Input         : zz    (Coefficients in zig-zag order)
%
%  Output        : coeffs(DCT coefficients in original order)
zigz =  [    1     2     6     7    15    16    28    29;
     3     5     8    14    17    27    30    43;
     4     9    13    18    26    31    42    44;
    10    12    19    25    32    41    45    54;
    11    20    24    33    40    46    53    55;
    21    23    34    39    47    52    56    61;
    22    35    38    48    51    57    60    62;
    36    37    49    50    58    59    63    64];
for k = 1:3
    coeffs_t = zz(zigz(:),k);
    coeffs_t = reshape(coeffs_t,[8,8]);
    coeffs(:,:,k) = coeffs_t;
end
end
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
function mv = searchmv(search,block)
best_ssd = inf;
    for i = 1:9
        for j = 1:9
            search_block = search(i:i+7,j:j+7);
            ssd_current = sum((search_block-block).^2,'all');
            if ssd_current<best_ssd
                best_ssd = ssd_current;
                mv = [i,j]-[5,5];
            end
        end
    end
end
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