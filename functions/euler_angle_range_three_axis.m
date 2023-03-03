function a = euler_angle_range_three_axis(angles)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% EULER_ANGLE_RANGE_THREE_AXIS Limits Euler angles range.
% For three-axis rotation, the angle ranges are [-pi, pi], [-pi/2, pi/2] and [-pi, pi]
% For two-axis rotation, the angle ranges are [-pi, pi], [0, pi] and [-pi, pi]!
%
% INPUT:
%   * angles,           Euler angles                                            (3 x 1) vector      [rad]
%
% OUTPUT:
%   * a,                Limited Euler angles                                    (3 x 1) vector      [rad]
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);

    if (~isequal(size(angles), [3 1]))
        error('angles must be a (3 x 1) vector.');
    end
    
    % convert the second angle in range [-pi, pi]
    a1 = angles(1);
    a2 = wrapToPi(angles(2));
    a3 = angles(3);
    %a1 = wrapToPi(a1);
    %a3 = wrapToPi(a3);
    
    % if the second angle is not within [-pi/2, pi/2]
    if a2 > pi / 2
        a1 = a1 + pi;
        a2 = pi - a2;
        a3 = a3 + pi;
    elseif a2 < - pi / 2
        a1 = a1 + pi;
        a2 = pi - a2;
        a3 = a3 + pi;
    end
    a1 = wrapToPi(a1);
    a3 = wrapToPi(a3);

    a = [a1; a2; a3];
    
end