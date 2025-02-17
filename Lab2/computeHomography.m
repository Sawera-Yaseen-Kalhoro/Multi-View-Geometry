function H12 = computeHomography(CL1uv,CL2uv, Model)
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

    switch (Model)

        case 'Translation'

            % rotation part
            R = eye(2);

            % translation part
            t_all = CL1uv-CL2uv;
            t = mean(t_all);
            t = t';

            % construct homography
            H12 = [R,t;
                0,0,1];
            
        case 'Similarity'

            % expressing the problem as Ax = B as in the lect slide 23
            % construct A and B
            for i = 1:length(CL1uv)
                A(2*i-1,:) = [CL2uv(i,1) -CL2uv(i,2) 1 0];
                A(2*i,:) = [CL2uv(i,2) CL2uv(i,1) 0 1];
                
                B(2*i-1,1) = CL1uv(i,1);
                B(2*i,1) = CL1uv(i,2);
            end

            x = (A'*A)\(A'*B);

            % % convert a and b into s and theta
            % theta = atan(x(2)/x(1));
            % s = x(1)/cos(theta);
            % 
            % % construct the homography
            % H12 = [s*cos(theta) -s*sin(theta) x(3);
            %     s*sin(theta) s*cos(theta) x(4);
            %     0 0 1];

            % construct the homography
            H12 = [x(1) -x(2) x(3);
                x(2) x(1) x(4);
                0 0 1];

            
        case 'Affine'

            % expressing the problem as Ax = B
            % construct A and B
            for i = 1:length(CL1uv)
                A(2*i-1,:) = [CL2uv(i,1) CL2uv(i,2) 1 0 0 0];
                A(2*i,:) = [0 0 0 CL2uv(i,1) CL2uv(i,2) 1];
                
                B(2*i-1,1) = CL1uv(i,1);
                B(2*i,1) = CL1uv(i,2);
            end

            x = (A'*A)\(A'*B);

            % construct the homography
            H12 = [x(1) x(2) x(3);
                x(4) x(5) x(6);
                0 0 1];

        
        case 'Projective'

            % expressing the problem as Ax = B
            % construct A and B
            for i = 1:length(CL1uv)
                A(2*i-1,:) = [CL2uv(i,1) CL2uv(i,2) 1 0 0 0 -CL2uv(i,1)*CL1uv(i,1) -CL2uv(i,2)*CL1uv(i,1)];
                A(2*i,:) = [0 0 0 CL2uv(i,1) CL2uv(i,2) 1 -CL2uv(i,1)*CL1uv(i,2) -CL2uv(i,2)*CL1uv(i,2)];
                
                B(2*i-1,1) = CL1uv(i,1);
                B(2*i,1) = CL1uv(i,2);
            end

            x = (A'*A)\(A'*B);

            % construct the homography
            H12 = [x(1) x(2) x(3);
                x(4) x(5) x(6);
                x(7) x(8) 1];

        
        otherwise
            warning('Invalid model, returning identity homography');
            H12 = eye(3);
            
    end
    
    
end





