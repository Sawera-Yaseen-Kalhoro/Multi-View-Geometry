function [F, sum_abs_diff, cam1_fig, cam2_fig] = compute_fundamental_matrix_normalization(x1,x2,P1,P2,F21,rank2)
% function to perform step 8-10 for normalized coordinates
% x1 = 2D points in image 1
% x2 = 2D points in image 2
% P1 = projection matrix for camera 1
% P2 = projection matrix for camera 2
% F21 = fundamental matrix obtained from camera parameters in step 4 
% (will be used as the ground truth)
% rank2 = bool value; enforces rank2 condition if true

imageSize = [256 256];
tx = -1000; %mm
ty = 190; %mm
tz = 230; %mm

% ensure that the points are in cartesian coordinate (not homogenous)
if size(x1,1) == 3
    x1 = [x1(1,:)./x1(3,:);x2(2,:)./x1(3,:)];
elseif size(x1,1) ~= 2
    fprintf("Invalid data points for x1.")
    return
end

if size(x2,1) == 3
    x2 = [x2(1,:)./x2(3,:);x2(2,:)./x2(3,:)];
elseif size(x2,1) ~= 2
    fprintf("Invalid data points for x2.")
    return
end

%% Normalization
% Compute the mean and standard deviation of x and y coordinates for both images
mean_x1 = mean(x1(1, :));
mean_y1 = mean(x1(2, :));
std_x1 = std(x1(1, :));
std_y1 = std(x1(2, :));

mean_x2 = mean(x2(1, :));
mean_y2 = mean(x2(2, :));
std_x2 = std(x2(1, :));
std_y2 = std(x2(2, :));

% Create transformation matrices T1 and T2 for coordinate normalization
T1 = [1/std_x1 0 -mean_x1/std_x1; 0 1/std_y1 -mean_y1/std_y1; 0 0 1];
T2 = [1/std_x2 0 -mean_x2/std_x2; 0 1/std_y2 -mean_y2/std_y2; 0 0 1];

% Apply coordinate normalization
x1_n = T1 * [x1; ones(1, size(x1, 2))];
x2_n = T2 * [x2; ones(1, size(x2, 2))];

%% step 8
% Compute the Fundamental Matrix Using 8- point method

U_n = [x2_n(1,:)'.*x1_n(1,:)' x2_n(1,:)'.*x1_n(2,:)' x2_n(1,:)' ...
 x2_n(2,:)'.*x1_n(1,:)' x2_n(2,:)'.*x1_n(2,:)' x2_n(2,:)' ...
 x1_n(1,:)' x1_n(2,:)' ones(length(x1_n),1) ];

fprintf("Condition number of U_n = %f\n", cond(U_n))

[U, D, V] = svd(U_n);

% Extract fundamental matrix from the column of V
% corresponding to the smallest singular value.
F = reshape(V(:,9),3,3)';

if rank2
    % Enforce rank2 constraint
    [U,D,V] = svd(F);
    F = U*diag([D(1,1) D(2,2) 0])*V';
end

%% aditional step for normalized coordinates
% invert F
F = T2'*F*T1;

%% step 9
% Compare the F21 (Step 4) and F (Step 8) matrices
% Normalize both matrices first
F21 = F21/F21(3,3);
F = F/F(3,3);

% Sum of Absolute Difference
sum_abs_diff = sum(abs(F21 - F),"all");

%% step 10
% Epipolar Geometry
% Draw 2D projections on image planes
cam1_fig = mvg_show_projected_points(x1(1:2,:),imageSize,'Image 1');
cam2_fig = mvg_show_projected_points(x2(1:2,:),imageSize,'Image 2');

% Draw epipolar lines
[~,~,c1_l_coeff,c2_l_coeff] = mvg_compute_epipolar_geom_modif(x1,x2,F);
[cam1_fig,cam2_fig] = mvg_show_epipolar_lines(cam1_fig, cam2_fig, c1_l_coeff,c2_l_coeff, [-400,1;300,400],'b');

% Compute epipoles
[~, ~, V] = svd(F');
ep1 = V(:, end);
ep1_normalized = ep1 / ep1(end);  % Normalize to make it homogeneous

[~, ~, U] = svd(F);
ep2 = U(:, end);
ep2_normalized = ep2 / ep2(end);  % Normalize to make it homogeneous

% Draw epipoles
[~,~] = mvg_show_epipoles(cam1_fig, cam2_fig,ep2_normalized,ep1_normalized);

end