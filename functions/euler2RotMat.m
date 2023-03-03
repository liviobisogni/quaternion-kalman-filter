function C = euler2RotMat(r)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% EULER2ROTMAT Converts Euler angles to rotation matrix (aka direction cosine matrix, DCM).
% The Euler angles rotate the frame n (navigation) to the frame b (body) according to 'zyx' sequence.
%
% INPUT:
%   * angles,           Euler angles (angles = [roll; pitch; yaw])              (3 x 1) vector      [rad]
%
% OUTPUT:
%   * C,                Coordinate transformation matrix from n to b            (3 x 3) matrix      []
%                       NB: That is v_b  = C * v_n.
%                           ('_b' or '_n' mean the vector 'v' is expressed in the frame b or n)
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(r), [3 1]))
        error('r must be a (3 x 1) vector.');
    end
    
    r1 = r(1);      % roll
    r2 = r(2);      % pitch
    r3 = r(3);      % yaw
    
%     r1 = -r1;
%     r2 = -r2;
%     r3 = -r3;
    
    C = rotX(r1) * rotY(r2) * rotZ(r3);                 % "old"; semi-OK
    
%     C = rotZ(r3) * rotY(r2) * rotX(r1);               % "new"; semi-OK!!!!!!
%     
%     C_n_b = rotX(r1) * rotY(r2) * rotZ(r3);
%     C_b_n = rotZ(r3) * rotY(r2) * rotX(r1);
%     
%     C_n_b_TRANS = C_n_b';
%     C_b_n_TRANS = C_b_n';
    

%     C = C';
    
%     C = rotZ(r1) * rotY(r2) * rotX(r3);                 % NEW, maybe
    
    
%         C = rotX(-r1) * rotY(-r2) * rotZ(-r3);         % CANC
    
%     C = rotX(r3) * rotY(r2) * rotZ(r1);         % CANC






%% ALTERNATIVE method
%     % [Eq. 1 "Quaternion-Based Unscented Kalman Filter for Accurate Indoor
%     %         Heading Estimation Using Wearable Multi-Sensor System" - Xuebing Yuan]
%     
%     cx = cos(-r1);       % cos(phi)
%     sx = sin(-r1);       % sin(phi)
%     cy = cos(-r2);       % cos(theta)
%     sy = sin(-r2);       % sin(theta)
%     cz = cos(-r3);       % cos(psi)
%     sz = sin(-r3);       % sin(psi)
%     
%     C = [cz * cy,                       sz * cy,                    -sy;
%          -sz * cx + cz * sy * sx,       cz * cx + sz * sy * sx,     cy * cx;
%          sz * sx + cz * cy * cx,        -cz * sx + sz * sy * cx,    cy * cx];
% 
%     C = C / norm(C);
    
end