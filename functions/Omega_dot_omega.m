function Omega__dot_omega = Omega_dot_omega(omega_next, omega_prev, dt)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 126 Trawny]
% OMEGA_DOT_OMEGA Defines the derivative of the turn rate (omega_dot) and the associated matrix
% (Omega(omega_dot)), which - in the linear case - is constant.
%
% INPUT:
%   * omega_next,       Angular velocity omega(k+1)                             (3 x 1) vector      [rad / s]
%   * omega_prev,       Angular velocity omega(k)                               (3 x 1) vector      [rad / s]
%   * dt,               Time step (dt = t(k+1) - t(k))                          scalar              [s]
%
% OUTPUT:
%   * Omega__dot_omega, Omega matrix                                            (3 x 3) matrix
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(3,3);

    if (~isequal(size(omega_next), [3 1]))
        error('omega_next must be a (3 x 1) vector.');
    end
    if (~isequal(size(omega_prev), [3 1]))
        error('omega_prev must be a (3 x 1) vector.');
    end
    if (~isscalar(dt))
        error('dt must be a scalar.');
    end
    
%     global dt

    % [Eq. 126 Trawny]
    Omega__dot_omega = Omega((omega_next - omega_prev) / dt);
    
end