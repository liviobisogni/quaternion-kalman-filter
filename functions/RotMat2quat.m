function q = RotMat2quat(C)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 98a - 98b - 99a - 99b Trawny]
% ROTMAT2QUAT Converts direction cosine matrix to quaternion.
%
% INPUT:
%   * C,                Coordinate transformation matrix                        (3 x 3) matrix      []
%
% OUTPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(C), [3 3]))
        error('C must be a (3 x 3) matrix.');
    end
%     if (norm(C) ~= 1)
%         warning('C should have norm 1.');
%     end
    if (abs(1 - norm(C)) > 1e-15)
        warning(['C should have norm 1:    norm(C) = ', num2str(norm(C)), ',    norm(C) - 1 = ', num2str(norm(C)-1)])
    end
        
    for i=1:3
        for j=1:3
            if (C(i,j) < -1)
                warning('C[%i,%i]=%f is smaller than -1', i, j, C(i,j));
            elseif (C(i,j) > 1)
                warning('C[%i,%i]=%f is greater than 1', i, j, C(i,j));
            end
        end
    end
    
    c11 = C(1,1);
    c12 = C(1,2);
    c13 = C(1,3);
    c21 = C(2,1);
    c22 = C(2,2);
    c23 = C(2,3);
    c31 = C(3,1);
    c32 = C(3,2);
    c33 = C(3,3);
    
    % [Eq. 96 Trawny]
    T = trace(C);
    
    % "Each of these four solutions will be singular when the pivotal
    % element is zero, but at least one will not be (since otherwise q
    % could not have unit norm)1. For maximum numerical accuracy, the form
    % with the largest pivotal element should be used."
    [~, max_index] = max([c11, c22, c33, T]);
    
    switch(max_index)
        
        case 1
            % [Eq. 98a Trawny]
            q1 = sqrt(1 + 2 * c11 - T) / 2;
            q = [q1;
                 (c12 + c21) / (4 * q1);
                 (c13 + c31) / (4 * q1);
                 (c23 - c32) / (4 * q1)];
            
        case 2
            % [Eq. 98b Trawny]
            q2 = sqrt(1 + 2 * c22 - T) / 2;
            q = [(c12 + c21) / (4 * q2);
                 q2;
                 (c23 + c32) / (4 * q2);
                 (c31 - c13) / (4 * q2)];
             
        case 3
            % [Eq. 99a Trawny]
            q3 = sqrt(1 + 2 * c33 - T) / 2;
            q = [(c13 + c31) / (4 * q3);
                 (c23 + c32) / (4 * q3);
                 q3;
                 (c12 - c21) / (4 * q3)];
             
        case 4
            % [Eq. 99b Trawny]
            q4 = sqrt(1 + T) / 2;
            q = [(c23 - c32) / (4 * q4);
                 (c31 - c13) / (4 * q4);
                 (c12 - c21) / (4 * q4);
                 q4];
             
        otherwise
            error('Unexpected behavior.')
    end
    
    % quaternion normalization, *** no need if dcm is really a dcm
    q = q / norm(q);

	
    % Ensure quaternion scalar part is non-negative
    if (q(4) < 0)
        q = -q;
    end
    
    q = [q(4); q(1); q(2); q(3)]; % poiché il nuovo paper usa lo scalare come primo elemento (anziché come ultimo)
    

    %% ALTERNATIVE method
% %     https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2015/01/matrix-to-quat.pdf
% %     Check also: % https://gunjanpatel.wordpress.com/quaternions-to-rotation-angles-how-to-use-them/
% 
%     if (c33 < 0)
%         if (c11 > c22)
%             t = 1 + c11 - c22 - c33;
%             q = [t; c12 + c21; c31 + c13; c23 - c32];
%         else
%             t = 1 - c11 + c22 - c33;
%             q = [c12 + c21; t; c23 + c32; c31 - c13];
%         end
%     else
%         if (c11 < -c22)
%             t = 1 - c11 - c22 + c33;
%             q = [c31 + c13; c23 + c32; t; c12 - c21];
%         else
%             t = 1 + m00 + m11 + m22;
%             q =[c23 - c32, c31 - c13, c12 - c21, t];
%         end
%     end
%     
%     q = q * 0.5 / sqrt(t);
%     
%     % ensure q(4) is non-negative
%     if (q(4) < 0)
%         q = -q;
%     end
%     
%     q = [q(4); q(1); q(2); q(3)]; % poiché il nuovo paper usa lo scalare come primo elemento (anziché come ultimo)

end