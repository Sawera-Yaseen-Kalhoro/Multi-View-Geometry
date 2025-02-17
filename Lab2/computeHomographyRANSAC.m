function H12 = computeHomographyRANSAC(CL1uv,CL2uv, Model)
%% computeHomography : estimate the Homography between two images according to Model 
%           Cluv1    set of points on image1 . Each row represents a 2-D point
%                       (u,v). Size: Nx2, with N number of points.
%           Cluv2    set of points on image2 . Each row represents a 2-D point
%                       (u,v). Size: Nx2, with N number of points.
%           Model       type of Homography to estimate. It has to be egual
%                       to one of the following strings: 'Translation',
%                       'Rigid', 'Similarity', 'Affine', 'Projective'.
%   Output
%           H12           estimated Homography of Model type. 3x3 matrix.
%

    % set some parameters for the RANSAC algorithm
    distance_threshold = 0.003;
    consensus_threshold = 0.8*length(CL1uv);

    switch (Model)

        case 'Translation'
            min_samples = 1;
            max_iter = 5;

        case 'Similarity'
            min_samples = 2;
            max_iter = 5;

        case 'Affine'
            min_samples = 3;
            max_iter = 7;

        case 'Projective'
            min_samples = 4;
            max_iter = 9;
    end

    % initialize variable to store the max # of consensus
    % in case the algorithm never reaches the consensus threshold
    max_consensus = 0;

    % initialize new lists for CL1uv and CL2uv containing only the inliers
    CL1uv_new = [];
    CL2uv_new = [];

    % apply RANSAC
    for iter = 1:max_iter
        % randomly select min number of samples and create a model
        rand_index = randi(length(CL1uv),min_samples,1);
        points1 = CL1uv(rand_index,:);
        points2 = CL2uv(rand_index,:);
        test_model = computeHomography(points1,points2,Model);

        % calculate consensus and store inliers
        projerror = projectionerrorvec(test_model,CL1uv,CL2uv)

        % initialize subsets of CL1uv and CL2uv that are in consensus
        S1 = [];
        S2 = [];

        for point_index = 1:length(CL1uv)
            if projerror(point_index) < distance_threshold
                S1 = [S1;CL1uv(point_index,:)];
                S2 = [S2;CL2uv(point_index,:)];
            end
        end
        num_consensus = length(S1);

        % if the number of consensus is greater than a threshold, 
        % re-create the model using the inliers and terminate
        if num_consensus > consensus_threshold
            CL1uv_new = S1;
            CL2uv_new = S2;
            break
        
        elseif num_consensus > max_consensus
            % if the number of consensus never reach the threshold, 
            % we need to store the biggest consensus
            max_consensus = num_consensus;
            S1_max = S1;
            S2_max = S2;
        end     
    end

    % if the number of consensus never reach the threshold,
    % use the points with the highest consensus
    if isempty(CL1uv_new)
        CL1uv_new = S1_max;
        CL2uv_new = S2_max;
    end

    % compute homography depending on the model chosen
    H12 = computeHomography(CL1uv_new,CL2uv_new, Model);

end