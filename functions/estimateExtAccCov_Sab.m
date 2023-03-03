function Q_hat_a_b = estimateExtAccCov_Sab(y_a)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. 37 Suh]
% ESTIMATEEXTACCACOV_SAB Accelerometer norm-based adaptive algorithm by A. M. Sabatini for estimating
% external acceleration covariance matrix Q__a_b.
%
% INPUT:
%   * y_a,              Measured acceleration                               (3 x 1) vector      [m/s^2]
%
% OUTPUT:
%   * Q_hat_a_b,        Estimated external acceleration covariance          (3 x 3) matrix      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
    if (~isequal(size(y_a), [3 1]))
        error('y_a must be a (3 x 1) vector.');
    end
    
    global g_standard
    global s                                                    % e.g., 1 or 10 or 100
    global epsilon_a                                            % e.g., 0.2                                 [m / s^2]
    global acc_norm_adapt_COUNTER                               % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!

    % [Eq. 37 Suh]
    if (abs(norm(y_a) - g_standard) < epsilon_a)                % No external acceleration detected
        Q_hat_a_b = zeros(3);
    else                                                        % External acceleration detected
        acc_norm_adapt_COUNTER = acc_norm_adapt_COUNTER + 1;    % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
        Q_hat_a_b = s * eye(3);
        
    end

end