function C_z = rotZ(alpha)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% ROTZ Basic rotation about z-axis by an angle alpha.
%
% INPUT:
%   * alpha,            angle                                                   scalar              [rad]
%
% OUTPUT:
%   * C,                Coordinate transformation matrix                        (3 x 3) matrix      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isscalar(alpha))
        error('alpha must be a scalar.');
    end
    
    ca = cos(alpha);
    sa = sin(alpha);
    
    C_z = [ca,  sa, 0;
           -sa, ca, 0;
           0,   0,  1];

end