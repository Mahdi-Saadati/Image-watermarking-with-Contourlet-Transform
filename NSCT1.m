clc;
clear all;
close all;

% Parameteters:
nlevels =[2,3];        % Decomposition level

Scaling_Factor=0.05;

InputImage=imread('Airplane.jpg');
InputImage=im2double(InputImage);
InputImage_R=InputImage(:,:,1);
InputImage_G=InputImage(:,:,2);
InputImage_B=InputImage(:,:,3);

watermark=imread('Peugeot_logo.jpg');
watermark=im2double(watermark);

% Nonsubsampled Contourlet decomposition
coeffs_R = nsctdec(InputImage_R,nlevels);
coeffs_G = nsctdec(InputImage_G,nlevels);
coeffs_B = nsctdec(InputImage_B,nlevels);

%disp(nlevels); disp(dfilter);

for i=1:4
      Energy_R(i)=(sum(sum(coeffs_R{1,2}{1,i})));
      Energy_G(i)=(sum(sum(coeffs_G{1,2}{1,i})));
      Energy_B(i)=(sum(sum(coeffs_B{1,2}{1,i})));
end

[MaxEnergy_R,j]=max(Energy_R);  %MaxEnergyy meghdar Maximum Energy zirtasavir ast va j shomare an zirtasvir ast
[MaxEnergy_G,k]=max(Energy_G);
[MaxEnergy_B,l]=max(Energy_B);

LL_R=coeffs_R{1,2}{1,j};
LL_G=coeffs_G{1,2}{1,k};
LL_B=coeffs_B{1,2}{1,l};

[LL_R_U,LL_R_S,LL_R_V]=svd(LL_R);
[LL_G_U,LL_G_S,LL_G_V]=svd(LL_G);
[LL_B_U,LL_B_S,LL_B_V]=svd(LL_B);

[m,n,~]=size(InputImage);
watermark=imresize(watermark,[m,n]);

watermark_R=watermark(:,:,1);
watermark_G=watermark(:,:,2);
watermark_B=watermark(:,:,3);

LL_R_S2=watermark_R*Scaling_Factor+LL_R_S;
LL_G_S2=watermark_G*Scaling_Factor+LL_G_S;
LL_B_S2=watermark_B*Scaling_Factor+LL_B_S;

LL_R_new=LL_R_U*LL_R_S2*LL_R_V';
LL_G_new=LL_G_U*LL_G_S2*LL_G_V';
LL_B_new=LL_B_U*LL_B_S2*LL_B_V';

coeffs_R{1,2}{1,j}=LL_R_new;
coeffs_G{1,2}{1,k}=LL_G_new;
coeffs_B{1,2}{1,l}=LL_B_new;

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
PSNR_value_Old=PSNR(InputImage,WatermarkingImage);
SSIM_value_Old=ssim(InputImage,WatermarkingImage);

%% Attacks image watermarking
%% 1) For Geometrical attack: cropping and rotation, you can use imcrop() and imrotate():

%WatermarkingImage = imrotate(WatermarkingImage,25,'bilinear','crop');    %Rotation(RO)
%WatermarkingImage = imcrop(WatermarkingImage);     %Crop (CR)
%WatermarkingImage=imresize(WatermarkingImage,0.10);  % Scaling (SC)
% WatermarkingImage=flipdim2(WatermarkingImage);      %Flips the columns, making a mirror image  %flip (FL).
% WatermarkingImage=flipdim1(WatermarkingImage);  %Flips the rows, making an upside-down image
%% 2) For Noising attack: Gaussian noise, you can use imnoise()

% WatermarkingImage = imnoise(WatermarkingImage,'gaussian');  %Gaussian Noise(GN)
% WatermarkingImage = imnoise(WatermarkingImage,'salt & pepper',0.05);   %Pepper and Salt noise (SP)
% WatermarkingImage = imnoise(WatermarkingImage,'speckle',0.05);  %Speckle Noise (SN)


%% 3) For Denoising attack: average filtering, you can use imfilter() or conv2()

% WatermarkingImage = medfilt3(WatermarkingImage);  %Median Filter (MF)

%% 4) For Format-compression attack: JPEG compression, you can use imwrite()

%  imwrite(WatermarkingImage,'C:\Compress_Image.jpeg','Mode','lossy','Quality',75);   %JPEG Compression (JPEG)
%  WatermarkingImage = imread('C:\Compress_Image.jpeg');
%  WatermarkingImage=im2double(WatermarkingImage);

%% 5) For Image-processing attack: histogram equalization (HE),contrast adjustment (CA), and gamma correction (GC), you can use histeq(), adapthisteq(), imadjust(), and intlut(), repectively.
% WatermarkingImage=imadjust(WatermarkingImage,0.001);   %Gamma Correction (GC)
% WatermarkingImage=histeq(WatermarkingImage);    %Histogram Equalization (HE)
% WatermarkingImage=imfilter_AverageFilter(WatermarkingImage,1);   %Average Filter (AF)
%WatermarkingImage=fspecial_GaussianFilter(WatermarkingImage,1);  %Gaussian low-pass Filter (GP)
%WatermarkingImage=fspecial_MutionFilter(WatermarkingImage,15,20);  %Motion Blur(MB)
% WatermarkingImage=imsharpen(WatermarkingImage);   %Sharpening (SH)
%%%watermarkingImage=skewing(WatermarkingImage,0.2,0.2);  %Shearing (SE)
% WatermarkingImage=fspecial_BlurringFilter(WatermarkingImage,0.2);  %Blurring (BL)


%% Color Image Watermark Extraction

WatermarkingImage_R=WatermarkingImage(:,:,1);
WatermarkingImage_G=WatermarkingImage(:,:,2);
WatermarkingImage_B=WatermarkingImage(:,:,3);

coeffs_R = nsctdec(WatermarkingImage_R,nlevels);
coeffs_G = nsctdec(WatermarkingImage_G,nlevels);
coeffs_B = nsctdec(WatermarkingImage_B,nlevels);

LL_R_new=coeffs_R{1,2}{1,j};
LL_G_new=coeffs_G{1,2}{1,k};
LL_B_new=coeffs_B{1,2}{1,l};

[LL_R_U,LL_R_S,LL_R_V]=svd(LL_R_new);
[LL_G_U,LL_G_S,LL_G_V]=svd(LL_G_new);
[LL_B_U,LL_B_S,LL_B_V]=svd(LL_B_new);

LL_R=LL_R_U*LL_R_S*LL_R_V';
LL_G=LL_G_U*LL_G_S*LL_G_V';
LL_B=LL_B_U*LL_B_S*LL_B_V';

coeffs_R{1,2}{1,j}=LL_R;
coeffs_G{1,2}{1,k}=LL_G;
coeffs_B{1,2}{1,l}=LL_B;

Extracting_Orginal_Image_R = nsctrec(coeffs_R);
Extracting_Orginal_Image_G = nsctrec(coeffs_G);
Extracting_Orginal_Image_B = nsctrec(coeffs_B);

Extracting_Orginal_Image=cat(3,Extracting_Orginal_Image_R,Extracting_Orginal_Image_G,Extracting_Orginal_Image_B);

ExtractWatermark_R=(LL_R_S2-LL_R_S)/Scaling_Factor;
ExtractWatermark_G=(LL_G_S2-LL_G_S)/Scaling_Factor;
ExtractWatermark_B=(LL_B_S2-LL_B_S)/Scaling_Factor;

ExtractWatermark=cat(3,ExtractWatermark_R,ExtractWatermark_G,ExtractWatermark_B);

% figure;
% subplot(1,3,1);
% imshow(watermark);
% title('Orginal watermark');
% subplot(1,3,2);
% imshow(ExtractWatermark);
% title('Extracting watermark');
% subplot(1,3,3);
% imshow(Extracting_Orginal_Image);
% title('Extracting Orginal Image');

PSNR_value_New=PSNR(watermark,ExtractWatermark);
SSIM_value_New=ssim(watermark,ExtractWatermark);
NC_value=(sum(sum(sum(watermark(:,:,:).*ExtractWatermark(:,:,:)))) / ((sqrt(sum(sum(sum(watermark(:,:,:).^2))))).*(sqrt(sum(sum(sum(ExtractWatermark(:,:,:).^2)))))));
CRC_value=(sum(sum(sum(watermark(:,:,:).*ExtractWatermark(:,:,:)))) / (sum(sum(sum(watermark(:,:,:).^2)))));

fprintf('PSNR=%.4f \nSSIM=%.4f \n\nPSNR=%.4f \nSSIM=%.4f  \nNC=%.4f  \nCRC=%.4f \n',PSNR_value_Old,SSIM_value_Old,PSNR_value_New,SSIM_value_New,NC_value,CRC_value)