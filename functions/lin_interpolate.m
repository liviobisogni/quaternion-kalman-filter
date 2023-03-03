function values_interpolated = lin_interpolate(T, values, N_samples, N_samples__theoretical)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. pre-34 Suh]
% LIN_INTERPOLATE Linearly interpolates values.
%
% INPUT:
%   * T,                        Time stamp array                            (1 x N_samples) vector  [s]
%   * values,                   Array with elements to be interpolated      (1 x N_samples) vector  [-]
%   * N_samples,                Number of samples of T and values           scalar                  []
%   * N_samples__theoretical,   Theoretical number of samples               scalar                  []
%   * dt,                       Constant time step                          scalar                  [s]
%
% OUTPUT:
%   * values_interpolated,    	N_samples__theoretical interpolated values  (1 x N_samples__theoretical) array  [-]
%
% N.B.: N_samples <= N_samples__theoretical
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    global dt
    
    % Check number of arguments
    narginchk(4,4);
    
    if (~isscalar(N_samples))
        error('N_samples must be a scalar.');
    end
    if (~isscalar(N_samples__theoretical))
        error('N_samples__theoretical must be a scalar.');
    end
    if (~isequal(size(T), [1 N_samples]))
        error('T must be a (1 x N_samples) vector.');
    end
    if (~isequal(size(values), [1 N_samples]))
        error('values must be a (1 x N_samples) vector.');
    end
    if (N_samples > N_samples__theoretical)
        error('N_samples__theoretical must be greater than N_samples.');
    end
    
%     dt = 0.1;                                                 % CANC
%     T       = [0.09; 0.21; 0.3; 0.4; 0.5; 0.6; 0.7; 0.9; 1];  % CANC
%     values  = [2;    5;    1;   1;   1;   1;   3;   1;   2];  % CANC
%     N_samples__theoretical = 10;                              % CANC

    T_theoretical = zeros(1, N_samples__theoretical);
    values_interpolated = zeros(1, N_samples__theoretical);

for j=1:(N_samples__theoretical - 0)
%     j                                   % DEBUG PRINT
    T_theoretical(j) = j * dt;          % {0.02, 0.04, 0.06, ...}
%     T_theoretical(j)                    % DEBUG PRINT
    
    % https://it.mathworks.com/matlabcentral/answers/578296-how-to-determine-the-two-closest-values-to-a-threshold-value
    % nearest value considering it is bigger or smaller than threshold:
    [T_lower, idx_lower] = max(T(T <= T_theoretical(j)));     	% nearest value before threshold
    if (abs(T_lower - T_theoretical(j)) < 1e-10)
%         disp('IF')
        values_interpolated(j) = values(idx_lower);
    else
        T_upper = min(T(T >= T_theoretical(j)));                % nearest value after threshold
        idx_upper = find(T == T_upper);

        % Linear interpolation
        values_interpolated(j) = values(idx_lower) + (T_theoretical(j) - T_lower) * (values(idx_upper) - values(idx_lower)) / (T_upper - T_lower);
%         values_interpolated(j)              % DEBUG PRINT

    end
end