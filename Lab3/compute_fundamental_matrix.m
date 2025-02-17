function [F, sum_abs_diff, cam1_fig, cam2_fig] = compute_fundamental_matrix(x1,x2,P1,P2,F21,rank2)
% function to perform step 8-10
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

%% step 8
% Compute the Fundamental Matrix Using 8- point method

U_n = [x2(1,:)'.*x1(1,:)' x2(1,:)'.*x1(2,:)' x2(1,:)' ...
 x2(2,:)'.*x1(1,:)' x2(2,:)'.*x1(2,:)' x2(2,:)' ...
 x1(1,:)' x1(2,:)' ones(length(x1),1) ];

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