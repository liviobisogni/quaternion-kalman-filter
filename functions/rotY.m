function C_y = rotY(alpha)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% ROTY Basic rotation about y-axis by an angle alpha.
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
    
    C_y = [ca,  0,  -sa;
           0,   1,  0;
           sa,  0,  ca];

end