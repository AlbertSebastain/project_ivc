
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
%im_encode = blockproc(im_quant,[4,4],@(block_sruct) CAVLC_encode(bock_struct.data));
%dst = im_encode(:);
im_zigz = blockproc(im_quant,[8,8],@(block_struct) ZigZag8x8(block_struct.data));
im_zigz = im_zigz(:);
dst = ZeroRunEnc_EoB(im_zigz,EOB);
end