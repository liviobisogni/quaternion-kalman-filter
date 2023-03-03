function Phi = Phi_matrix(y_g__next, Phi_mode)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. post-15 Suh] ('precise') OR [Eq. 16a Suh] ('approximated') OR [Eq. 187 Trawny] ('Trawny')
% PHI_MATRIX Computes the transition matrix Phi(t + dt, t).
%
% INPUT:
%   * y_g__next,        Measured angular velocity (i.e., gyro output) y_g(k)    (3 x 1) vector      [rad / s]
%   * dt,               Time step (dt = t(k+1) - t(k))                          scalar              [s]
%   * Phi_mode,         Possible values:    * 'precise':        [Eq. post-15 Suh]
%                                           * 'approximated':   [Eq. 16a Suh]
%                                           * 'Trawny':         [Eq. 187 Trawny]
%
% OUTPUT:
%   * Phi,              Transition matrix Phi(t + dt, t)                        (9 x 9) matrix
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(2,2);
    
    global dt
    
    if (~isequal(size(y_g__next), [3 1]))
        error('y_g__next must be a (3 x 1) vector.');
    end
%     if (~isscalar(dt))
%         error('dt must be a scalar.');
%     end
    if (~ischar(Phi_mode))
        error('Phi_mode must be a string.');
    end
    
    % [Eq. post-9 Suh]
    A = [-skewSymmetric(y_g__next),     - 1/2 * eye(3),     zeros(3);
         zeros(3),                      zeros(3),           zeros(3);
         zeros(3),                      zeros(3),           zeros(3)];
    
    if strcmpi(Phi_mode, 'Trawny')              % "old" Phi_mode
        % [Eq. 187 Trawny]
        Phi = [Theta(y_g__next, dt),    Psi(y_g__next, dt);
               zeros(3),                eye(3)];
    elseif strcmpi(Phi_mode, 'precise')
        % [Eq. post-15 Suh]
        Phi = exp(A * dt);
    elseif strcmpi(Phi_mode, 'approximated')
        % [Eq. 16a Suh]
        Phi = eye(9) + A * dt + 1/2 * A^2 * dt^2;
    else
        error('Invalid mode for Phi (Phi_mode). Possible values: ''precise'', ''approximated'', and ''Trawny''.');
    end
    
end