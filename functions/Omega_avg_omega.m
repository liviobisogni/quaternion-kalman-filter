function Omega__avg_omega = Omega_avg_omega(omega_next, omega_prev, dt)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 129 Trawny]
% OMEGA_AVG_OMEGA Computes the Omega(omega_avg) matrix.
%
% INPUT:
%   * omega_next,       Angular velocity omega(k+1)                             (3 x 1) vector      [rad / s]
%   * omega_prev,       Angular velocity omega(k)                               (3 x 1) vector      [rad / s]
%   * dt,               Time step (dt = t(k+1) - t(k))                          scalar              [s]
%
% OUTPUT:
%   * Omega__avg_omega, Omega(omega_avg) matrix                                 (3 x 3) matrix      [rad / s]
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

    % [Eq. 129 Trawny]
    Omega__avg_omega = Omega(omega_next) + 1/2 * Omega_dot_omega(omega_next, omega_prev, dt) * dt;
    
end