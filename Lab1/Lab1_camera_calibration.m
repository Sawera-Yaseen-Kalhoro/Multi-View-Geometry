clc; clear all; close all

%% step 1
% Intrinsic parameters
au = 557.0943; av = 712.9824;
u0 = 326.3819; v0 = 298.6679;
% Location of the world reference frame in camera coordinates in mm
Tx = 100; Ty = 0; Tz = 1500;
% World rotation w.r.t. camera coordinates
% Euler XYX1 angles
Phix = 0.8*pi/2;
Phiy = -1.8*pi/2;
Phix1 = pi/5;

%% step 2
% intrinsic matrix
K = [au 0 u0;
    0 av v0;
    0 0 1];

% rotation matrices
cRw_x = [1 0 0;
    0 cos(Phix) -sin(Phix);
    0 sin(Phix) cos(Phix)];

cRw_y = [cos(Phiy) 0 sin(Phiy);
    0 1 0;
    -sin(Phiy) 0 cos(Phiy)];

cRw_x1 = [1 0 0;
    0 cos(Phix1) -sin(Phix1);
    0 sin(Phix1) cos(Phix1)];

cRw = cRw_x*cRw_y*cRw_x1;

% translation matrix
ctw = [Tx;Ty;Tz];

% homogenous transformation matrix
cTw = [cRw ctw;
    0 0 0 1];

% projection matrix
P = K*[eye(3) zeros(3,1)]*cTw;

% normalize for later comparisons
P_norm = P/P(3,4);

%% step 3
% define 6 points in 3D
p1 = [-17;374;-400];
p2 = [-171;413;390];
p3 = [213;-152;216];
p4 = [-17;-446;329];
p5 = [344;6;-130];
p6 = [-8;-53;80];

% gather all the points
p = [p1,p2,p3,p4,p5,p6];

% alt: generate 6 new points
%p = generate_random_3d_points(6);

% plot for checking purposes
plot3(p(1,:),p(2,:),p(3,:), 'ro')
grid on

%% step 4
% compute projection
for i = 1:6
    pp_hom(:,i) = P*[p(:,i);1];
end

%% step 5
% convert from homogenous to cartesian
pp_cart = [pp_hom(1,:)./pp_hom(3,:);
    pp_hom(2,:)./pp_hom(3,:)];

% plot the 2D points
figure
plot(pp_cart(1,:),pp_cart(2,:), 'ro')
grid on

%% step 6
% estimate the camera parameters

% create the Q matrix and B vector
for i = 1:6
    Q(2*i-1,:) = [p(1,i) p(2,i) p(3,i) 1 0 0 0 0 -pp_cart(1,i)*p(1,i) -pp_cart(1,i)*p(2,i) -pp_cart(1,i)*p(3,i)];
    Q(2*i,:) = [0 0 0 0 p(1,i) p(2,i) p(3,i) 1 -pp_cart(2,i)*p(1,i) -pp_cart(2,i)*p(2,i) -pp_cart(2,i)*p(3,i)];

    B(2*i-1,1) = pp_cart(1,i);
    B(2*i,1) = pp_cart(2,i);
end

% compute A vector
A = Q\B;

% create the new projection matrix using the elements of A
P_new = [A(1) A(2) A(3) A(4);
    A(5) A(6) A(7) A(8);
    A(9) A(10) A(11) 1];

%% step 7
% extract intrinsic parameters
[K_new,cRw_new] = get_intrinsics_from_proj_matrix(P_new);

%% step 8
% add noise
noise = normrnd(0,0.5,size(pp_cart));
pp_cart_noisy = pp_cart + noise;

% repeat step 6
% create the Q matrix and B vector
for i = 1:6
    Q_noisy(2*i-1,:) = [p(1,i) p(2,i) p(3,i) 1 0 0 0 0 -pp_cart_noisy(1,i)*p(1,i) -pp_cart_noisy(1,i)*p(2,i) -pp_cart_noisy(1,i)*p(3,i)];
    Q_noisy(2*i,:) = [0 0 0 0 p(1,i) p(2,i) p(3,i) 1 -pp_cart_noisy(2,i)*p(1,i) -pp_cart_noisy(2,i)*p(2,i) -pp_cart_noisy(2,i)*p(3,i)];

    B_noisy(2*i-1,1) = pp_cart_noisy(1,i);
    B_noisy(2*i,1) = pp_cart_noisy(2,i);
end

% compute A vector
A_noisy = Q_noisy\B_noisy;

% create the new projection matrix using the elements of A
P_new_noisy = [A_noisy(1) A_noisy(2) A_noisy(3) A_noisy(4);
    A_noisy(5) A_noisy(6) A_noisy(7) A_noisy(8);
    A_noisy(9) A_noisy(10) A_noisy(11) 1];

% extract intrinsic parameters
[K_noisy,cRw_noisy] = get_intrinsics_from_proj_matrix(P_new_noisy);

%% step 9
% compute projection
for i = 1:6
    pp_noisy_hom(:,i) = P_new_noisy*[p(:,i);1];
end

% convert into cartesian
pp_noisy_cart = [pp_noisy_hom(1,:)./pp_noisy_hom(3,:);
    pp_noisy_hom(2,:)./pp_noisy_hom(3,:)];

% compute projection error
proj_error = 0;
for i = 1:6
    proj_error = proj_error + sqrt((pp_cart(1,i)-pp_noisy_cart(1,i))^2 + (pp_cart(2,i)-pp_noisy_cart(2,i))^2);
end
proj_error = proj_error/6;

%% step 10
% For simplicity, we have compiled step 8-9 into a function called
% avg_projection_error. Thus, we are just calling that function here.

% compute for 10 points
p_10 = [p,generate_random_3d_points(4)];
proj_error(2) = avg_projection_error(P,p_10)

% compute for 50 points
p_50 = [p_10,generate_random_3d_points(40)];
proj_error(3) = avg_projection_error(P,p_50)

% plot the average projection errors
number_of_points = [6,10,50];
figure
plot(number_of_points,proj_error, 'o-')
grid on