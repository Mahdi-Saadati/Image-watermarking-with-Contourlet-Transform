function f29(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)

imwrite(WatermarkingImage,'K:\Compress_Image.jpeg','Mode','lossy','Quality',80);   %JPEG Compression (JPEG)
WatermarkingImage = imread('K:\Compress_Image.jpeg');     %JPEG Compression (JPEG)
WatermarkingImage=im2double(WatermarkingImage); 

Watermark_Extraction(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B);

end