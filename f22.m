function f22(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)

WatermarkingImage=fspecial_MutionFilter(WatermarkingImage,15,45);

Watermark_Extraction(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B);

end