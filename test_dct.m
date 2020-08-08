img = imread('foreman0020.bmp');
img = double(img);
img_rgb = img;
img = ictRGB2YCbCr(img);
qScale = 1;
%im_dct = blockproc(img,[4,4],@(block_struct) DCT4X4(block_struct.data));
%im_quant = blockproc(im_dct,[4,4],@(block_struct) Quant4x4(block_struct.data,qScale));
%image_dequant = blockproc(im_quant,[4,4],@(block_struct) DeQuant4x4(block_struct.data,qScale));
%im_idct = blockproc(image_dequant,[4,4],@(block_struct) IDCT4x4(block_struct.data));
%im_dct = blockproc(img,[4,4],@(block_struct) DCT4X4(block_struct.data));
 im_dct = blockproc(img,[8,8],@(block_struct) DCT8x8(block_struct.data));
 im_quant = blockproc(im_dct,[8,8],@(block_struct) Quant8x8(block_struct.data,qScale));
 image_dequant = blockproc(im_quant,[8,8],@(block_struct) DeQuant8x8(block_struct.data,qScale));
 im_idct = blockproc(image_dequant,[8,8],@(block_struct) IDCT8x8(block_struct.data));
 best_a = 0;
 best_b = 0;
 best_PSNR = 0;
 for a = 1:20;
     for b = 1:20;
         
        deblock_im = deblock(im_idct,a,b);
        im_rgb_deblock = ictYCbCr2RGB(deblock_im);
        PSNR_deb = calcPSNR(im_rgb_deblock,img_rgb);
        if best_PSNR < PSNR_deb
            best_PSNR = PSNR_deb;
            best_a = a;
            best_b = b;
        end
     end
 end
im_idct = ictYCbCr2RGB(im_idct);
PSNR = calcPSNR(im_idct,img_rgb);
%PSNR_deb = calcPSNR(im_rgb_deblock,img_rgb);
im_idct = uint8(im_idct);
imshow(im_idct);