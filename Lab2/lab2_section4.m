%% main execution file for lab 2 section 4 - Estimating a Homography
clc, clear all, close all

% add path to all folders
addpath("MiamiSet00","matlab_helper_functions_lab2","DataSet01")

%% step 1 -> created a new function computeHomography.m

%% step 2
load("Features.mat")

% homography from image 00 to 01
image1filename = "00.png";
image2filename = "01.png";
fprintf('Homography for images: %s - %s',image1filename,image2filename)

CL1uv = Features(1).xy;
CL2uv = Features(2).xy;

I1 = imread(image1filename);
I2 = imread(image2filename);

Model = 'Similarity'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Translation'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Projective'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Affine'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

% homography from image 00 to 02
image2filename = "02.png";
fprintf('Homography for images: %s - %s',image1filename,image2filename)

CL2uv = Features(3).xy;

I2 = imread(image2filename);

Model = 'Similarity'
H12 = computeHomography(CL1uv,CL2uv, Model)
% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Translation'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Projective'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Affine'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

% homography from image 00 to 03
image2filename = "03.png";
fprintf('Homography for images: %s - %s',image1filename,image2filename)

CL2uv = Features(4).xy;

I2 = imread(image2filename);

Model = 'Similarity'
H12 = computeHomography(CL1uv,CL2uv, Model)
% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Translation'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Projective'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

Model = 'Affine'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

%% step 3
% Miami coral reef
image1filename = 'imgl01311.jpg';
image2filename = 'imgl01396.jpg';

% generate paired points using matchsiftmodif.m
distRatio = 0.5;
drawMatches = false;
[CL1uv,CL2uv] = matchsiftmodif(image1filename, image2filename, distRatio, drawMatches);

% compute homography
fprintf('Homography for images: %s - %s',image1filename,image2filename)
I1 = imread(image1filename);
I2 = imread(image2filename);

Model = 'Similarity'
H12 = computeHomography(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);