function a_b = create_extAcc(a_b_OLD, t_i, ext_acc, length)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% CREATE_EXTACC Generates constant external acceleration ext_acc from t_i to t_i + length, overwriting
% previous external acceleration a_b_OLD.
%
% INPUT:
%   * a_b_OLD,          Previous external acceleration; it'll be overwritten    (3 x N_samples) vector  [m / s^2]
%   * t_i,              Ext. acc. starting (initial) time                       scalar                  [s]
%   * ext_acc,          Ext. acc. constant value                                (3 x 1) vector          [m / s^2]
%   * length,           Ext. acc. duration                                      scalar                  [s]
%
% OUTPUT:
%   * a_b,              External acceleration a_b(t)                            (3 x N_samples) vector  [m / s^2]
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    global dt N_samples
%     global f_sampling                                               % [Hz]
    f_sampling = 1 / dt;                                               % [Hz]

    % Check number of arguments
    narginchk(4,4);
    
    if (~isequal(size(a_b_OLD), [3 N_samples]))
        error('a_b_OLD must be a (3 x N_samples) vector.');
    end
    if (~isscalar(t_i))
        error('t_i must be a scalar.');
    end
    if (~isequal(size(ext_acc), [3 1]))
        error('ext_acc must be a (3 x 1) vector.');
    end
    if (~isscalar(length))
        error('length must be a scalar.');
    end

    % a_b = zeros(3, N_samples);
    a_b = a_b_OLD;

    for i = t_i * f_sampling : t_i * f_sampling + length * f_sampling
        a_b(:,i) = ext_acc;
    end
    
end