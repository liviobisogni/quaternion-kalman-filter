function C_from_q = quat2RotMat(q)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 90 Trawny] ([Eq. 1 Suh])
% QUAT2ROTMAT Computes the rotation matrix associated to the quaternion q.
%
% INPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      []
%
% OUTPUT:
%   * C_from_q,         Coordinate transformation matrix associated to q        (3 x 3) matrix      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(q), [4 1]))
        error('q must be a (4 x 1) vector.');
    end
%     if (norm(q) ~= 1)
%         warning('q should be a unit vector.');
%     end
    if (abs(1 - norm(q)) > 1e-13)   % 1e-15
        warning(['q should be a unit vector:    norm(q) = ', num2str(norm(q)), ',    norm(q) - 1 = ', num2str(norm(q)-1)])
    end
    
    q = q / norm(q);                                                         % CANC?
    
    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);
    
    %% [Eq. 79 Trawny] (OLD)
%     C_from_q = Csi_R(q)' * Psi_L(q);


    %% [Eq. 1 Suh] (ALTERNATIVE formula; it DOES work)
%
%     C_from_q = [2 * q0^2 + 2 * q1^2 - 1,     2 * q1 * q2 + 2 * q0 * q3,      2 * q1 * q3 - 2 * q0 * q2;
%                 2 * q1 * q2 - 2 * q0 * q3,   2 * q0^2 + 2 * q2^2 - 1,        2 * q2 * q3 + 2 * q0 * q1;
%                 2 * q1 * q3 + 2 * q0 * q2,   2 * q2 * q3 - 2 * q0 * q1,      2 * q0^2 + 2 * q3^2 - 1];
            

    %% [Eq. 90 Trawny] (apparently more precise)

    q4 = q0;            % poiché io / Suh usiamo il primo elemento come scalare di un quaternione, mentre Trawny il quarto
    q1q1 = q1 * q1;
    q1q2 = q1 * q2;
    q1q3 = q1 * q3;
    q1q4 = q1 * q4;
    q2q2 = q2 * q2;
    q2q3 = q2 * q3;
    q2q4 = q2 * q4;
    q3q3 = q3 * q3;
    q3q4 = q3 * q4;
    q4q4 = q4 * q4;
    
    C_from_q = [q1q1 - q2q2 - q3q3 + q4q4,   2 * (q1q2 + q3q4),              2 * (q1q3 - q2q4);
                2 * (q1q2 - q3q4),           - q1q1 + q2q2 - q3q3 + q4q4,    2 * (q2q3 + q1q4);
                2 * (q1q3 + q2q4),           2 * (q2q3 - q1q4),              - q1q1 - q2q2 + q3q3 + q4q4];

            
    % C_from_q = C_from_q / norm(C_from_q);     % CANC?

end