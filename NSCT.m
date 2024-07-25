tic
clc;
clear all;
close all;

% Parameteters:
Scaling_Factor=0.05;

watermark0=imread('Numbers.jpg');
watermark0=im2double(watermark0);

InputImage=imread('Tree.jpg');
InputImage=im2double(InputImage);
[m1,n1,~]=size(InputImage);

InputImage_R=InputImage(:,:,1);
InputImage_G=InputImage(:,:,2);
InputImage_B=InputImage(:,:,3);

watermark=imread('Numbers.jpg');
watermark=im2double(watermark);
[m2,n2,~]=size(watermark);

watermark=imresize(watermark,[m1,n1]);

watermark_R=watermark(:,:,1);
watermark_G=watermark(:,:,2);
watermark_B=watermark(:,:,3);

% Nonsubsampled Contourlet decomposition
nlevels =[2,3];        % Decomposition level

coeffs_R = nsctdec(InputImage_R,nlevels);
coeffs_G = nsctdec(InputImage_G,nlevels);
coeffs_B = nsctdec(InputImage_B,nlevels);

LL_R=coeffs_R{1,1};
LL_G=coeffs_G{1,1};
LL_B=coeffs_B{1,1};

[LL_R_U,LL_R_S,LL_R_V]=svd(LL_R);
[LL_G_U,LL_G_S,LL_G_V]=svd(LL_G);
[LL_B_U,LL_B_S,LL_B_V]=svd(LL_B);

LL_R_S2=watermark_R*Scaling_Factor+LL_R_S;
LL_G_S2=watermark_G*Scaling_Factor+LL_G_S;
LL_B_S2=watermark_B*Scaling_Factor+LL_B_S;

LL_R_new=LL_R_U*LL_R_S2*LL_R_V';
LL_G_new=LL_G_U*LL_G_S2*LL_G_V';
LL_B_new=LL_B_U*LL_B_S2*LL_B_V';

coeffs_R{1,1}=LL_R_new;
coeffs_G{1,1}=LL_G_new;
coeffs_B{1,1}=LL_B_new;

WatermarkingImage_R = nsctrec(coeffs_R);
WatermarkingImage_G = nsctrec(coeffs_G);
WatermarkingImage_B = nsctrec(coeffs_B);

WatermarkingImage=cat(3,WatermarkingImage_R,WatermarkingImage_G,WatermarkingImage_B);

% figure;
% subplot(1,2,1);
% imshow(InputImage);
% title('Orginal Image');
% subplot(1,2,2);
% imshow(WatermarkingImage);
% title('Watermarking Image');

%% PSNR & SSIM
% PSNR_value_Old=PSNR(InputImage,WatermarkingImage);
% SSIM_value_Old=ssim(InputImage,WatermarkingImage);

%% Attacks image watermarking
f0(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Non Attack
f1(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Average Filter 3*3(AF)
f2(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Average Filter 5*5(AF)
f3(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Gaussian low-pass Filter 3*3(GP)
f4(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Gaussian low-pass Filter 5*5(GP)
f5(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Median Filter 3*3(MF)
f6(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Median Filter 5*5(MF)
f7(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %%Crop [150,150](CR)
f8(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Rotation 5(RO)
f9(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B) %Rotation 30(RO)
f10(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Rotation 45(RO)
f11(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Rotation 70(RO)
f12(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Rotation 110(RO)
f13(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Translation [5,10](TR)
f14(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Translation [10,10](TR)
f15(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Translation [10,15](TR)
f16(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Flip horizontal(FL)
f17(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Flip vertical(FL)
f18(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Blurring 0.3(BL)
f19(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Blurring 0.5(BL)
f20(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Blurring 1(BL)
f21(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Motion Blur (15,20)(MB)
f22(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Motion Blur (15,45)(MB)
f23(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Sharpening 0.3(SH)
f24(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Sharpening 0.5(SH)
f25(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Sharpening 1(SH)
f26(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%JPEG Compression 5(JPEG)
f27(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%JPEG Compression 10(JPEG)
f28(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%JPEG Compression 20(JPEG)
f29(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%JPEG Compression 80(JPEG)
f30(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%JPEG Compression 90(JPEG)
f31(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gaussian Noise 0.001(GN)
f32(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gaussian Noise 0.01(GN)
f33(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gaussian Noise 0.1(GN)
f34(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gaussian Noise 0.3(GN)
f35(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gamma Correction 0.001(GC)
f36(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gamma Correction 0.01(GC)
f37(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gamma Correction 0.1(GC)
f38(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Gamma Correction 0.3(GC)
f39(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Pepper and Salt noise 0.001(SP)
f40(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Pepper and Salt noise 0.01(SP)
f41(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Pepper and Salt noise 0.1(SP)
f42(WatermarkingImage,watermark0,LL_R_S2,LL_G_S2,LL_B_S2,Scaling_Factor,m2,n2,nlevels,LL_R,LL_G,LL_B)%Pepper and Salt noise 0.3(SP)

%% Color Image Watermark Extraction

% WatermarkingImage_R=WatermarkingImage(:,:,1);
% WatermarkingImage_G=WatermarkingImage(:,:,2);
% WatermarkingImage_B=WatermarkingImage(:,:,3);
% 
% coeffs_R = nsctdec(WatermarkingImage_R,nlevels);
% coeffs_G = nsctdec(WatermarkingImage_G,nlevels);
% coeffs_B = nsctdec(WatermarkingImage_B,nlevels);
% 
% LL_R_new=coeffs_R{1,1};
% LL_G_new=coeffs_G{1,1};
% LL_B_new=coeffs_B{1,1};
% 
% [LL_R_U,LL_R_S,LL_R_V]=svd(LL_R_new);
% [LL_G_U,LL_G_S,LL_G_V]=svd(LL_G_new);
% [LL_B_U,LL_B_S,LL_B_V]=svd(LL_B_new);
% 
% % LL_R=LL_R_U*LL_R_S*LL_R_V';
% % LL_G=LL_G_U*LL_G_S*LL_G_V';
% % LL_B=LL_B_U*LL_B_S*LL_B_V';
% 
% coeffs_R{1,1}=LL_R;
% coeffs_G{1,1}=LL_G;
% coeffs_B{1,1}=LL_B;
% 
% Extracting_Orginal_Image_R = nsctrec(coeffs_R);
% Extracting_Orginal_Image_G = nsctrec(coeffs_G);
% Extracting_Orginal_Image_B = nsctrec(coeffs_B);
% 
% Extracting_Orginal_Image=cat(3,Extracting_Orginal_Image_R,Extracting_Orginal_Image_G,Extracting_Orginal_Image_B);
% 
% ExtractWatermark_R=(LL_R_S2-LL_R_S)/Scaling_Factor;
% ExtractWatermark_G=(LL_G_S2-LL_G_S)/Scaling_Factor;
% ExtractWatermark_B=(LL_B_S2-LL_B_S)/Scaling_Factor;
% 
% ExtractWatermark=cat(3,ExtractWatermark_R,ExtractWatermark_G,ExtractWatermark_B);
% ExtractWatermark=imresize(ExtractWatermark,[m2,n2]);
% 
% PSNR_value_New=psnr(watermark0,ExtractWatermark);
% SSIM_value_New=ssim(watermark0,ExtractWatermark);
% NC_value =(sum(sum(sum(watermark0(:,:,:).*ExtractWatermark(:,:,:)))) / ((sqrt(sum(sum(sum(watermark0(:,:,:).^2))))).*(sqrt(sum(sum(sum(ExtractWatermark(:,:,:).^2)))))));
% fprintf('%.4f \n%.4f  \n%.4f \n',PSNR_value_New,SSIM_value_New,NC_value)

% figure;
% subplot(1,3,1);
% imshow(watermark0);
% title('Orginal watermark');
% subplot(1,3,2);
% imshow(ExtractWatermark);
% title('Extracting watermark');
% subplot(1,3,3);
% imshow(Extracting_Orginal_Image);
% title('Extracting Orginal Image');

toc
