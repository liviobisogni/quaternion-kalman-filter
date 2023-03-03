function q = euler2quat(angles)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [https://marc-b-reynolds.github.io/math/2017/04/18/TaitEuler.html#mjx-eqn%3Aeq%3Atait,
%  Eq. 2 - 3 - 4 - 5]
% EULER2QUAT Converts Euler angles to quaternion.
%
% INPUT:
%   * angles,           Euler angles (angles = [roll; pitch; yaw])              (3 x 1) vector      [rad]
%
% OUTPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(angles), [3 1]))
        error('angles must be a (3 x 1) vector.');
    end
    
%     angles = wrapToPi(angles);
    
    c1 = cos(angles(1) / 2);
    c2 = cos(angles(2) / 2);
    c3 = cos(angles(3) / 2);
    s1 = sin(angles(1) / 2);
    s2 = sin(angles(2) / 2);
    s3 = sin(angles(3) / 2);
    
    
    % [https://marc-b-reynolds.github.io/math/2017/04/18/TaitEuler.html#mjx-eqn%3Aeq%3Atait, Eq. 2 - 3 - 4 - 5]
    q = [c1 * c2 * c3 + s1 * s2 * s3;
         s1 * c2 * c3 - c1 * s2 * s3;
         c1 * s2 * c3 + s1 * c2 * s3;
         c1 * c2 * s3 - s1 * s2 * c3];
        
    
% IT GIVES THE SAME RESULTS (YOU CAN DELETE IT):    
% %      void zyx_to_quat(quat_t* q, vec3_t* v)
% % {
%    hy = 0.5 * angles(3);     % v->z
%    hp = 0.5 * angles(2);     % v->y
%    hr = 0.5 * angles(1);     % v->x
%    ys = sin(hy);
%    yc = cos(hy);
%    ps = sin(hp);
%    pc = cos(hp);
%    rs = sin(hr);
%    rc = cos(hr);
% 
%   q = [rc * pc * yc + rs * ps * ys;
%        rs * pc * yc - rc * ps * ys;
%   	   rc * ps * yc + rs * pc * ys;
%   	   rc * pc * ys - rs * ps * yc];
% % }

     
    % Normalize the quaternion
    q = q / norm(q);
    
    % Ensure quaternion scalar part is non-negative
    if (q(1) < 0)
        q = -q;
    end
            
end