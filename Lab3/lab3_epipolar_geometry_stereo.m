%% Main execution script for Lab 3
clear all; clc; close all;
%addpath("mvg_functions")

%% step 1
% define camera 1
au1 = 100;
av1 = 120;
gamma = 0;
uo1 = 128;
vo1 = 128;
imageSize = [256 256];

%% step 2
% define camera 2
au2 = 90;
av2 = 110;
uo2 = 128;
vo2 = 128;

% transformation from world coordinates to camera 2
%xyz euler
ax = 0.1; by = pi/4; cz = 0.2;
%translation
tx = -1000; %mm
ty = 190; %mm
tz = 230; %mm

% define rotations and translations from world frame to each camera
% camera 1
wRc1 = eye(3);
wtc1 = [0 0 0]';

% camera 2
wRc2 = [1 0 0; 0 cos(ax) -sin(ax); 0 sin(ax) cos(ax)] * ...
    [cos(by) 0 sin(by);0 1 0; -sin(by) 0 cos(by)] * ...
    [cos(cz) -sin(cz) 0; sin(cz) cos(cz) 0; 0 0 1];
wtc2 = [tx ty tz]';

%% step 3
% compute camera intrinsics
K1 = [au1 gamma uo1;
    0 av1 vo1;
    0 0 1];

K2 = [au2 gamma uo2;
    0 av2 vo2;
    0 0 1];

% transformation from camera 1 to camera 2
c1Tw = inv([wRc1 wtc1; 0 0 0 1]);
wTc2 = [wRc2 wtc2; 0 0 0 1];
c1Tc2 = c1Tw*wTc2;

% projection matrices for both cameras
c2Tw = inv(wTc2);

P1 = K1*[eye(3) zeros(3,1)]*c1Tw;
P2 = K2*[eye(3) zeros(3,1)]*c2Tw;

%% step 4
% split c1Tc2
R12 = c1Tc2(1:3,1:3); % rotation matrix from camera 1 to 2
t12 = c1Tc2(1:3,4); % translation vector from camera 1 to 2

% construct skew-symmetric matrix
t12x = [0 -t12(3) t12(2);
    t12(3) 0 -t12(1);
    -t12(2) t12(1) 0];

F21 = inv(K2')*R12'*t12x*inv(K1);

%% step 5
% define object points
X(:,1) = [100;-400;2000];
X(:,2) = [300;-400;3000];
X(:,3) = [500;-400;4000];
X(:,4) = [700;-400;2000];
X(:,5) = [900;-400;3000];
X(:,6) = [100;-50;4000];
X(:,7) = [300;-50;2000];
X(:,8) = [500;-50;3000];
X(:,9) = [700;-50;4000];
X(:,10) = [900;-50;2000];
X(:,11) = [100;50;3000];
X(:,12) = [300;50;4000];
X(:,13) = [500;50;2000];
X(:,14) = [700;50;3000];
X(:,15) = [900;50;4000];
X(:,16) = [100;400;2000];
X(:,17) = [300;400;3000];
X(:,18) = [500;400;4000];
X(:,19) = [700;400;2000];
X(:,20) = [900;400;3000];

%% step 6
% compute projected images
x1 = P1*[X; ones(1,length(X))];
x1 = [x1(1,:)./x1(3,:); x1(2,:)./x1(3,:)]; % convert from hom to cart

x2 = P2*[X; ones(1,length(X))];
x2 = [x2(1,:)./x2(3,:); x2(2,:)./x2(3,:)]; % convert from hom to cart

%% step 7
% Using Helper Function
% Plotting the projected images without noise
cam1_fig_no_noise = mvg_show_projected_points(x1, imageSize, "Image 1 - No Noise");
cam2_fig_no_noise = mvg_show_projected_points(x2, imageSize, "Image 2 - No Noise");

%% step 8
% Compute the Fundamental Matrix Using 8- point method

U_n = [x2(1,:)'.*x1(1,:)' x2(1,:)'.*x1(2,:)' x2(1,:)' ...
 x2(2,:)'.*x1(1,:)' x2(2,:)'.*x1(2,:)' x2(2,:)' ...
 x1(1,:)' x1(2,:)' ones(length (X),1) ];

[U, D, V] = svd(U_n);

% Extract fundamental matrix from the column of V
% corresponding to the smallest singular value.
F = reshape(V(:,9),3,3)';

% Enforce rank2 constraint
[U,D,V] = svd(F);
F = U*diag([D(1,1) D(2,2) 0])*V';

%% step 9
% Compare the F21 (Step 4) anf F (Step 8) matrices
% Normalize both matrices first
F21 = F21/F21(3,3);
F = F/F(3,3);

% Sum of Absolute Difference
sum_abs_diff = sum(abs(F21 - F),"all");
fprintf("No noise\n")
fprintf('Sum of Absolute Difference between F21 and F: %f\n', sum_abs_diff);

%% step 10
% Epipolar Geometry
[~, ~, c1_l_coeff_no_noise, c2_l_coeff_no_noise] = mvg_compute_epipolar_geom_modif(x1, x2, F);
[cam1_fig_no_noise, cam2_fig_no_noise] = mvg_show_epipolar_lines(cam1_fig_no_noise, cam2_fig_no_noise, c1_l_coeff_no_noise, c2_l_coeff_no_noise, [-400, 1; 300, 400], 'b');

% Compute epipoles
[~, ~, V] = svd(F');
ep1 = V(:, end);
ep1_normalized = ep1 / ep1(end);  % Normalize to make it homogeneous

[~, ~, U] = svd(F);
ep2 = U(:, end);
ep2_normalized = ep2 / ep2(end);  % Normalize to make it homogeneous


% Draw epipoles
[~,~] = mvg_show_epipoles(cam1_fig_no_noise, cam2_fig_no_noise,ep2_normalized,ep1_normalized);

%% Step 11 
% Add noise
noise_x1 = normrnd (0, 0.5, size (x1));
noise_x2 = normrnd (0 ,0.5, size (x2));

x1_noise1 = x1 + noise_x1;
x2_noise1 = x2 + noise_x2;

%% Step 12
% Repeat steps 8-10 with noisy 2D points using the created function
fprintf("Coordinates with Noise [-1,+1]\n")
[F_noise1, sum_abs_diff_noise1, cam1_fig_noise1, cam2_fig_noise1] = compute_fundamental_matrix(x1_noise1,x2_noise1,P1,P2,F21,false);

fprintf('Sum of Absolute Difference between F21 and F_noise1: %f\n', sum_abs_diff_noise1);

% If we impose the rank2 condition
[F_noise1, sum_abs_diff_noise1, cam1_fig_noise1, cam2_fig_noise1] = compute_fundamental_matrix(x1_noise1,x2_noise1,P1,P2,F21,true);

%% Step 13
% Increase the noise
increased_noise_x1 = normrnd (0, 1, size (x1));
increased_noise_x2 = normrnd (0 ,1, size (x2));

x1_noise2 = x1 + increased_noise_x1;
x2_noise2 = x2 + increased_noise_x2;

% repeat steps 8-12 using the created function
fprintf("Coordinates with Noise [-2,+2]\n")
[F_noise2, sum_abs_diff_noise2, cam1_fig_noise2, cam2_fig_noise2] = compute_fundamental_matrix(x1_noise2,x2_noise2,P1,P2,F21,false);

fprintf('Sum of Absolute Difference between F21 and F_noise2: %f\n', sum_abs_diff_noise2);

% If we impose rank2 condition
[F_noise2, sum_abs_diff_noise2, cam1_fig_noise2, cam2_fig_noise2] = compute_fundamental_matrix(x1_noise2,x2_noise2,P1,P2,F21,true);

%% step 14
% coordinate data normalization on the noisy data
fprintf("Normalized Coordinates\n")

% normalize coordinates and repeat steps 8-12 using the created function
[F_normalized, sum_abs_diff_normalized, cam1_fig_normalized, cam2_fig_normalized] = compute_fundamental_matrix_normalization(x1_noise2,x2_noise2,P1,P2,F21,false);

fprintf('Sum of Absolute Difference between F21 and F_normalized: %f\n', sum_abs_diff_normalized);