function proj_error = avg_projection_error(P,points)
    % P = projection matrix to convert 3D into 2D points (before adding
    % noise)
    % points = a collection of 3D points, size = 3xn

    n = length(points); % number of points
    % compute projection
    for i = 1:n
        pp_hom(:,i) = P*[points(:,i);1];
    end
    
    % convert from homogenous to cartesian
    pp_cart = [pp_hom(1,:)./pp_hom(3,:);
        pp_hom(2,:)./pp_hom(3,:)];
    
    % % plot the 2D points
    % figure
    % plot(pp_cart(1,:),pp_cart(2,:), 'ro')
    % grid on
    
    % add noise
    noise = normrnd(0,0.5,size(pp_cart));
    pp_cart_noisy = pp_cart + noise;
    
    % create the Q matrix and B vector
    for i = 1:n
        Q_noisy(2*i-1,:) = [points(1,i) points(2,i) points(3,i) 1 0 0 0 0 -pp_cart_noisy(1,i)*points(1,i) -pp_cart_noisy(1,i)*points(2,i) -pp_cart_noisy(1,i)*points(3,i)];
        Q_noisy(2*i,:) = [0 0 0 0 points(1,i) points(2,i) points(3,i) 1 -pp_cart_noisy(2,i)*points(1,i) -pp_cart_noisy(2,i)*points(2,i) -pp_cart_noisy(2,i)*points(3,i)];
    
        B_noisy(2*i-1,1) = pp_cart_noisy(1,i);
        B_noisy(2*i,1) = pp_cart_noisy(2,i);
    end
    
    % compute A vector
    A_noisy = (Q_noisy'*Q_noisy)\(Q_noisy'*B_noisy);
    
    % create the new projection matrix using the elements of A
    P_new_noisy = [A_noisy(1) A_noisy(2) A_noisy(3) A_noisy(4);
        A_noisy(5) A_noisy(6) A_noisy(7) A_noisy(8);
        A_noisy(9) A_noisy(10) A_noisy(11) 1];
    
    % extract intrinsic parameters
    [K_noisy,cRw_noisy] = get_intrinsics_from_proj_matrix(P_new_noisy);
    
    % compute new (noisy) projection
    for i = 1:n
        pp_noisy_hom(:,i) = P_new_noisy*[points(:,i);1];
    end
    
    % convert into cartesian
    pp_noisy_cart = [pp_noisy_hom(1,:)./pp_noisy_hom(3,:);
        pp_noisy_hom(2,:)./pp_noisy_hom(3,:)];
    
    % compute projection error
    proj_error = 0;
    for i = 1:n
        proj_error = proj_error + sqrt((pp_cart(1,i)-pp_noisy_cart(1,i))^2 + (pp_cart(2,i)-pp_noisy_cart(2,i))^2);
    end
    proj_error = proj_error/n;
end