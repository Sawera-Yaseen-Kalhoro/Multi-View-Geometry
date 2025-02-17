% Helper script for lab 2 of MVG section 2
clear all;
close all;
clc;

addpath("MiamiSet00\","matlab_helper_functions_lab2\")

image1filename = 'imgl01311.jpg';
image2filename = 'imgl01396.jpg';

[image1, descriptors1, loc1] = siftprecomputed(image1filename);
showkeys(image1, loc1)

[image2, descriptors2, loc2] = siftprecomputed(image2filename);
showkeys(image2, loc2)


