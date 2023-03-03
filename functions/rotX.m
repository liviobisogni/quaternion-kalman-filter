function C_x = rotX(alpha)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% ROTX Basic rotation about x-axis by an angle alpha.
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
    
    C_x = [1,   0,      0;
           0,   ca,     sa;
           0,   -sa,    ca];

end