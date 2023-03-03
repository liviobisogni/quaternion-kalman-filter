function q_inverse = quatInverse(q)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 13 Trawny]
% QUATINVERSE Takes the inverse of a given quaternion
% (N.B.: The inverse rotation is described by the inverse - or complex conjugate for unit quaternions -
% quaternion)
%
% INPUT:
%   * q,                Quaternion q(k) (q = [q0; q1; q2; q3], q0 is the scalar)        (4 x 1) vector      []
%
% OUTPUT:
%   * q_inverse,        (not unitary) Quaternion, inverse of q                          (4 x 1) vector      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(q), [4 1]))
        error('q must be a (4 x 1) vector.');
    end

    q_inverse = [-q(1);
                 -q(2);
                 -q(3);
                  q(4)];
              
    q_inverse = q_inverse / norm(q_inverse)^2;

end