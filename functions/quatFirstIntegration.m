function q_next = quatFirstIntegration(q_prev, omega_next, omega_prev, integration_mode)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 18 Suh] ('Suh') OR [Eq. 131 Trawny] ('Trawny') OR [Eq. 8 Yuan] ('Yuan')
% QUATFIRSTINTEGRATION First order quaternion integrator.
% It makes the assumption of a linear evolution of omega during the integration interval dt.
% NB: DO use 'Suh'; DO NOT use 'Trawny' !!!
%
% INPUT:
%   * q,                Quaternion q(k) (q = [q0; q1; q2; q3], q0 is the scalar)(4 x 1) vector      []
%   * omega_next,       Angular velocity omega(k+1)                             (3 x 1) vector      [rad/s]
%   * omega_prev,       Angular velocity omega(k)                               (3 x 1) vector      [rad/s]
%   * dt,               Integration interval (dt = t(k+1) - t(k))               scalar              [s]
%   * integration_mode, Modality of integration:    * 'Suh',    [Eq. 18 Suh]
%                                                   * 'Trawny', [Eq. 131 Trawny]
%                                                   * 'Yuan',   [Eq. 8 Yuan]
%
% OUTPUT:
%   * q_next,           q(k+1); unit quaternion from integration of q           (4 x 1) vector      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(4,4);
    
    global dt
    
    if (~isequal(size(q_prev), [4 1]))
        error('q_prev must be a (4 x 1) vector.');
    end
    if (~isequal(size(omega_next), [3 1]))
        error('omega_next must be a (3 x 1) vector.');
    end
    if (~isequal(size(omega_prev), [3 1]))
        error('omega_prev must be a (3 x 1) vector.');
    end
%     if (~isscalar(dt))
%         error('dt must be a scalar.');
%     end
    if (~ischar(integration_mode))
        error('integration_mode must be a string or a char array. Possible values: ''Trawny'' or ''Suh'' or ''Yuan''');
    end

%     global dt

    if strcmpi(integration_mode, 'Trawny')      % Trawny
    %% OLD
    % [Eq. 131 Trawny]
    a = 0.5 * Omega_avg_omega(omega_next, omega_prev, dt) * dt;
    b = Omega(omega_next) * Omega(omega_prev) - Omega(omega_prev) * Omega(omega_next);
    q_next = (expm(a) + 1/48 * b * dt^2) * q_prev;

    elseif strcmpi(integration_mode, 'Suh')     % Suh
    %% NEW
    % [Eq. 18 Suh]
    q_next = (eye(4) + 3/4 * Omega(omega_next) * dt - 1/4 * Omega(omega_prev) * dt ...
              - 1/6 * norm(omega_next)^2 * dt^2 - 1/24 * Omega(omega_next) * Omega(omega_prev) * dt^2 ...
              - 1/48 * norm(omega_next)^2 * Omega(omega_next) * dt^3) * q_prev;

    elseif strcmpi(integration_mode, 'Yuan')    % Yuan
        % [Eq. 8 Yuan "Quaternion-Based Unscented Kalman Filter for Accurate Indoor
        %              Heading Estimation Using Wearable Multi-Sensor System" - Xuebing Yuan]
        w1 = omega_next(1);
        w2 = omega_next(2);
        w3 = omega_next(3);
        w1 = (omega_next(1) + omega_prev(1)) / 2;   % aggiunte mie
        w2 = (omega_next(2) + omega_prev(2)) / 2;   % aggiunte mie
        w3 = (omega_next(3) + omega_prev(3)) / 2;   % aggiunte mie
        A = 0.5 * Omega(omega_next);
        theta = sqrt((w1*dt)^2 + (w2 * dt)^2 + (w3 * dt)^2);
        if (theta ~= 0)
            q_next = (cos(theta / 2) * eye(4) + sin(theta / 2) / (theta / 2) * A * dt) * q_prev;
        else
            q_next = (cos(theta / 2) * eye(4) + 1 * A * dt) * q_prev;   % sin(x)/x ---> 1 (for x--->0)
        end
                  
    else
        error('Invalid integration mode. Possible values: ''Trawny'' or ''Suh'' or ''Yuan''');
    end

    %% IMPLEMENTAZIONE MENO PRECISA (CANC):
    % https://ahrs.readthedocs.io/en/latest/filters/angular.html
    % q_next = (eye(4) + 0.5 * Omega_avg_omega(omega_next, omega_prev) * dt) * q_next;

    %% Altra implementazione
%     q_next = quatMultiplication(q_prev, exp(0.5 * omega2quat(omega_next) * dt)); %% PROVALO
%     q(t0+dt) = q(t0)*exp((1/2)*W*dt);


    q_next = q_next / norm(q_next);   % CANC?

    % Ensure quaternion scalar part is non-negative
    if (q_next(1) < 0)
        q_next = -q_next;
    end

end