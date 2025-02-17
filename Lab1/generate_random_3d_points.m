function points = generate_random_3d_points(num_points)
    % Define the range for the coordinates
    min_range = -480;
    max_range = 480;

    % Generate and gather num_points random 3D points within the specified range
    points = (max_range - (min_range)) * rand(3, num_points) + min_range;
end
