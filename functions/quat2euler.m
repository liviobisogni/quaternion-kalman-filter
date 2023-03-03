function r = quat2euler(q)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [https://marc-b-reynolds.github.io/math/2017/04/18/TaitEuler.html#mjx-eqn%3Aeq%3Atait]
% QUAT2EULER Converts quaternion to Euler angles.
%
% INPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      []
%
% OUTPUT:
%   * r,                Euler angles (r = [roll; pitch; yaw])                   (3 x 1) vector      [rad]
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(q), [4 1]))
        error('q must be a (4 x 1) vector.');
    end
    
    % Normalize the quaternion
    q = q / norm(q);

%   x=q->x, y=q->y, z=q->z, w=q->w;

    q0 = q(1);
    q1 = q(2);
    q2 = q(3);
    q3 = q(4);

    t0 = q1 * q1 - q3 * q3;
    t1 = q0 * q0 - q2 * q2;
    xx = 0.5 * (t0 + t1);                               % 1/2 x of x'
    xy = q1 * q2 + q0 * q3;                             % 1/2 y of x'
    xz = q0 * q2 - q1 * q3;                             % 1/2 z of x'
    t  = xx * xx + xy * xy;                             % cos(theta)^2
    yz = 2 * (q2 * q3 + q0 * q1);                       % z of y'

    r3 = atan2(xy, xx);                                 % yaw   (psi)       v->z
    r2 = atan(xz/sqrt(t));                              % pitch (theta)     v->y

    if (t ~= 0)
        r1 = atan2(yz, t1 - t0);                        % rol   (phi)       v->x
    else
        r1 = (2 * atan2(q1,q0) - sgnd(xz) * r3);        % rol   (phi)       v->x
        disp('See QUAT2EULER.M!!!')
    end
    
    r = [r1; r2; r3];
    
%     r = euler_angle_range_three_axis(r);
    
end