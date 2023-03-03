function skew_symmetric_matrix = skewSymmetric(p)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. post-9 Suh]
% SKEW_SIMMETRIC computes the skew-symmetric matrix operator [qx], formed from p.
%
% INPUT:
%   * p,                  Possible values:    * angular velocity,               (3 x 1) vector      [rad / s]
%                                           * quaternion,                       (4 x 1) vector      []
%
% OUTPUT:
%	* skew_symmetric_matrix,      Skew-symmetric matrix formed from p           (3 x 3) matrix      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(p), [4 1]) && ~isequal(size(p), [3 1]))
        error('p must be a (4 x 1) or a (3x1) vector.');
    end

    if (size(p, 1) == 3)        % p is an angular velocity
        p1 = p(1);
        p2 = p(2);
        p3 = p(3);
    elseif (size(p, 1) == 4)    % p is a quaternion
        p1 = p(2);
        p2 = p(3);
        p3 = p(4);
    else
        error('Wrong array size.');
    end

    % [Eq. post-9 Suh]
    skew_symmetric_matrix = [0     -p3      p2;
                             p3     0      -p1;
                            -p2     p1      0];
                   
                                     
%     % DO NOT use this:
%     % [Eq. 14, "A New Quaternion-Based Kalman Filter for Real-Time Attitude
%                 Estimation Using the Two-Step Geometrically-Intuitive
%                 Correction Algorithm" - Kaiqiang Feng]
%
%     skew_symmetric_matrix = [0      p3     -p2;
%                              -p3     0       p1;
%                              p2     -p1      0];

end