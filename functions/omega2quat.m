function q = omega2quat(omega)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% OMEGA2QUAT Creates a pure imaginary quaternion (i.e., null scalar part, that is, q(1) = 0) from the
% angular velocity omega.
% Please note: the quaternion thus obtained is NOT a unitary quaternion.
%
% INPUT:
%   * omega,            Angular velocity omega(k)                               (3 x 1) vector      [rad / s]
%
% OUTPUT:
%   * q,                Quaternion (q = [q0; q1; q2; q3], q0 is the scalar)     (4 x 1) vector      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(omega), [3 1]))
        error('omega must be a (3 x 1) vector.');
    end
    
    q = [0; omega];
    
end
