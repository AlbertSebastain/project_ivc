function dst = IntraDecode(image, img_size , qScale,EOB)
%  Function Name : IntraDecode.m
%  Input         : image (zero-run encoded image, 1xN)
%                  img_size (original image size)
%                  qScale(quantization scale)
%  Output        : dst   (decoded image)
    image_rund = ZeroRunDec_EoB(image,EOB);
    image_rund = reshape(image_rund, [img_size(1)/4*16,img_size(2)/4*img_size(3)]);
    image_dezigz = blockproc(image_rund, [16,3],@(block_struct) DeZigZag4x4(block_struct.data));
    image_dequant = blockproc(image_dezigz,[4,4],@(block_struct) DeQuant4x4(block_struct.data,qScale));
    image_idc = blockproc(image_dequant,[4,4],@(block_struct) IDCT4x4(block_struct.data));
    %dst = ictYCbCr2RGB(image_idc);
    dst = image_idc;
end