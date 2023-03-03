function Q_d = Q_d_matrix(y_g__next, Q_d_mode)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. post-15 Suh] ('precise') OR [Eq. 16b Suh] ('approximated') OR [Eq. 208 Trawny] ('Trawny')
% Q_D_MATRIX Computes the Noise Covariance Matrix Q_d.
%
% INPUT:
%   * y_g__next,        Measured angular velocity (i.e., gyro output) y_g(k)    (3 x 1) vector      [rad / s]
%   * dt,               Time step (dt = t(k+1) - t(k))                          scalar              [s]
%   * time_old,         Current time stamp (t(k+1))                             scalar              [s]
%   * time_new,         Previous time stamp (t(k))                              scalar              [s]
%   * Q_d_mode,         Possible values:    * 'precise':        [Eq. post-15 Suh]
%                                           * 'approximated':   [Eq. 16a Suh]
%
% OUTPUT:
%   * Q_d,                Noise Covariance Matrix Q_d(k)                        (9 x 9) matrix
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

% function Q_d = Q_d_matrix(y_g__next, dt, time_old, time_new, Q_d_mode)

    % Check number of arguments
    narginchk(2,2);
    
    if (~isequal(size(y_g__next), [3 1]))
        error('y_g__next must be a (3 x 1) vector.');
    end
%     if (~isscalar(dt))
%         error('dt must be a scalar.');
%     end
%     if (~isscalar(time_old))
%         error('time_old must be a scalar.');
%     end
%     if (~isscalar(time_new))
%         error('time_new must be a scalar.');
%     end
    if (~ischar(Q_d_mode))
        error('Q_d_mode must be a string or a char array.');
    end
    
    global dt
%     global sigma_r sigma_w OMEGA_SMALL_THRESHOLD    % CANC ???????
    global Q
    
    A = [-skewSymmetric(y_g__next),     - 1/2 * eye(3),     zeros(3);
         zeros(3),                      zeros(3),           zeros(3);
         zeros(3),                      zeros(3),           zeros(3)];
    
%                                           * 'Trawny':         [Eq. 187 Trawny]
%     if strcmpi(Q_d_mode, 'Trawny')  % old mode; DO NOT use it! Btw it gives a 6x6 matrix, not a 9x9 matrix
%         
%         if (norm(y_g__next) < OMEGA_SMALL_THRESHOLD)
%             % [Eq. 212 Trawny]
%             Q11 = sigma_r^2 * eye(3) + sigma_w^2 * (eye(3) * dt^3 / 3 + ...
%                  2 * dt^5 / factorial(5) * skewSymmetricSquared(y_g__next));
%             % [Eq. 213 Trawny]
%             Q12 = -sigma_w^2 * (eye(3) * dt^2 / 2 - dt^3 / factorial(3) * skewSymmetric(y_g__next) + ...
%                   dt^4 / factorial(4) * skewSymmetricSquared(y_g__next));
%         else
%             % [Eq. 209 Trawny]
%             a = ((norm(y_g__next) * dt)^3 / 3 + 2 * sin(norm(y_g__next) * dt) - ...
%                 2 * norm(y_g__next) * dt) / norm(y_g__next)^5;
%             Q11 = sigma_r^2 * dt * eye(3) + sigma_w^2 * (eye(3) * dt^3 / 3 + ...
%                 a * skewSymmetricSquared(y_g__next));
%             % [Eq. 210 Trawny]
%             b = (norm(y_g__next) * dt - sin(norm(y_g__next) * dt)) / norm(y_g__next)^3;
%             c = ((norm(y_g__next) * dt)^2 / 2 + cos(norm(y_g__next) * dt) - 1) / norm(y_g__next)^4;
%             Q12 = -sigma_w^2 * (eye(3) * dt^2 / 2 + b * skewSymmetric(y_g__next) + c * skewSymmetricSquared(y_g__next));
%         end
%         % [Eq. 211 Trawny]
%         Q22 = sigma_w^2 * dt * eye(3);
%         % [Eq. 208 Trawny]
%         Q_d = [Q11,     Q12;
%                Q12',    Q22];
% 
%     elseif strcmpi(Q_d_mode, 'precise')

           
    if strcmpi(Q_d_mode, 'precise')
        
        % [Eq. post-15 Suh]
        fun = @(t) expm(A .* t) * Q * expm(A .* t)';
        Q_d = integral(fun, time_old, time_new, 'ArrayValued', true);
        
    elseif strcmpi(Q_d_mode, 'approximated')
        
        % [Eq. 16b Suh]
        Q_d = Q * dt + 1/2 * A * Q + 1/2 * Q * A';
                
    else
        error('Invalid mode for Q_d (Q_d_mode). Possible values: ''precise'' and ''approximated''.');
    end
       
end