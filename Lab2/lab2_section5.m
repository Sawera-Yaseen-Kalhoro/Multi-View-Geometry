%% main execution file for lab 2 section 5 - Improving the Registration Accuracy
clc, clear all, close all

% add path to all folders
addpath("MiamiSet00","matlab_helper_functions_lab2","DataSet01")

%% step 1 -> created a new function computeHomographyRANSAC.m

%% step 2
% Miami coral reef
image1filename = 'imgl01311.jpg';
image2filename = 'imgl01396.jpg';

% generate paired points using matchsiftmodif.m
distRatio = 0.6;
drawMatches = false;
[CL1uv,CL2uv] = matchsiftmodif(image1filename, image2filename, distRatio, drawMatches);

% compute homography
fprintf('Homography for images: %s - %s',image1filename,image2filename)
I1 = imread(image1filename);
I2 = imread(image2filename);

Model = 'Similarity'
H12 = computeHomographyRANSAC(CL1uv,CL2uv, Model)

% This is a fuction that warps image I2 into the frame of image I1 and
% shows the result with red and green colors
figure; showwarpedimages(I1,I2,H12);

% Compare with the result before RANSAC
H12_before = computeHomography(CL1uv,CL2uv, Model)
figure; showwarpedimages(I1,I2,H12_before);

% Compute reprojection error
errorVec = projectionerrorvec(H12,CL1uv,CL2uv)