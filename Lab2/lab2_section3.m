%% main execution file for lab 2 section 3 - Feature Association
clc, clear all, close all

% add path to all folders
addpath("MiamiSet00","matlab_helper_functions_lab2","DataSet01")

%% step 1

image1filename = 'imgl01311.jpg';
image2filename = 'imgl01396.jpg';

% Do feature association with the modified match.m function
distRatio = [0.4, 0.6, 0.8];
drawMatches = true;

for i = 1:3
    [CL1uv,CL2uv] = matchsiftmodif(image1filename, image2filename, distRatio(i), drawMatches);
end

%% step 2

load("H12.mat")
drawMatches = false; % suppress the figures
distThreshold = 50;

% Initialize stacks for plotting. 
% This can be implemented more elegantly using matlab's table structure
distRatioStack = [];
avgRepErrorStack = [];
maxRepErrorStack = [];
numPointsStack = [];
numPointsOverThresholdStack = []; 

for distRatio = 0.3:0.05:0.9
    [CL1uv,CL2uv] = matchsiftmodif(image1filename, image2filename, distRatio, drawMatches);
    fprintf('Projection error for distRatio = %.2f', distRatio)
    errorVec = projectionerrorvec(H12,CL1uv,CL2uv)

    if isempty(errorVec)
        disp('No features were associated');
    else
        distRatioStack = [distRatioStack; distRatio];
        avgRepErrorStack = [avgRepErrorStack; mean(errorVec)];
        maxRepErrorStack = [maxRepErrorStack; max(errorVec)];
        numPointsStack = [numPointsStack; size(CL1uv,1)];
        numPointsOverThresholdStack = [numPointsOverThresholdStack; sum(errorVec > distThreshold)];
    end
end

% Draw two versions of the same plot, one with linear y-axis and the other
% with log y-axis, to better see the small vs large values
figure
subplot(1,2,1)
plot(distRatioStack,avgRepErrorStack,distRatioStack,maxRepErrorStack,distRatioStack,numPointsStack,distRatioStack,numPointsOverThresholdStack)
xlabel('Distance Ratio')
ylabel('Error and Num Matches (linear scale)')
legend('Avg Rep Error','Max Rep Error','Num Matches','Num matches > dist threshold')

subplot(1,2,2)
semilogy(distRatioStack,avgRepErrorStack,distRatioStack,maxRepErrorStack,distRatioStack,numPointsStack,distRatioStack,numPointsOverThresholdStack)
xlabel('Distance Ratio')
ylabel('Error and Num Matches (log scale)')
legend('Avg Rep Error','Max Rep Error','Num Matches','Num matches > dist threshold')