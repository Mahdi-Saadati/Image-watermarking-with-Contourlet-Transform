function Watermark_Extraction(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)
global  ExtractWatermark;
% global  m2 n2; 
global  LL_R_S  LL_G_S  LL_B_S;
% global  LL_R_S2 LL_G_S2 LL_B_S2;

WatermarkingImage_R=WatermarkingImage(:,:,1);
WatermarkingImage_G=WatermarkingImage(:,:,2);
WatermarkingImage_B=WatermarkingImage(:,:,3);

coeffs_R = nsctdec(WatermarkingImage_R,nlevels);
coeffs_G = nsctdec(WatermarkingImage_G,nlevels);
coeffs_B = nsctdec(WatermarkingImage_B,nlevels);

LL_R_new=coeffs_R{1,1};
LL_G_new=coeffs_G{1,1};
LL_B_new=coeffs_B{1,1};

[LL_R_U,LL_R_S,LL_R_V]=svd(LL_R_new);
[LL_G_U,LL_G_S,LL_G_V]=svd(LL_G_new);
[LL_B_U,LL_B_S,LL_B_V]=svd(LL_B_new);

% LL_R=LL_R_U*LL_R_S*LL_R_V';
% LL_G=LL_G_U*LL_G_S*LL_G_V';
% LL_B=LL_B_U*LL_B_S*LL_B_V';

coeffs_R{1,1}=LL_R;
coeffs_G{1,1}=LL_G;
coeffs_B{1,1}=LL_B;

Extracting_Orginal_Image_R = nsctrec(coeffs_R);
Extracting_Orginal_Image_G = nsctrec(coeffs_G);
Extracting_Orginal_Image_B = nsctrec(coeffs_B);

Extracting_Orginal_Image=cat(3,Extracting_Orginal_Image_R,Extracting_Orginal_Image_G,Extracting_Orginal_Image_B);

ExtractWatermark_R=(LL_R_S2-LL_R_S)/Scaling_Factor;
ExtractWatermark_G=(LL_G_S2-LL_G_S)/Scaling_Factor;
ExtractWatermark_B=(LL_B_S2-LL_B_S)/Scaling_Factor;

ExtractWatermark=cat(3,ExtractWatermark_R,ExtractWatermark_G,ExtractWatermark_B);
ExtractWatermark=imresize(ExtractWatermark,[m2,n2]);


PSNR_value_New=PSNR(watermark0,ExtractWatermark);
SSIM_value_New=ssim(watermark0,ExtractWatermark);
NC_value =(sum(sum(sum(watermark0(:,:,:).*ExtractWatermark(:,:,:)))) / ((sqrt(sum(sum(sum(watermark0(:,:,:).^2))))).*(sqrt(sum(sum(sum(ExtractWatermark(:,:,:).^2)))))));
fprintf('%.4f \n%.4f  \n%.4f \n',PSNR_value_New,SSIM_value_New,NC_value)
% fprintf('%.4f\n',Best_pos)

end