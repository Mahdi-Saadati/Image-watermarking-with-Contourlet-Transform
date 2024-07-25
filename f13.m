function f13(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)

WatermarkingImage = imtranslate(WatermarkingImage,[5,10],'OutputView','same');

Watermark_Extraction(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B);
end