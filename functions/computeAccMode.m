function mode = computeAccMode(lambda, mu)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. pre-34 Suh]
% COMPUTEACCMODE Computes Shu's acceleration mode; it is able to detect external acceleration.
%
% INPUT:
%   * lambda,               The 3 eigenvalues of the matrix U, at times k, k-1, k-2     (3 x (M_2+1)) matrix
%   * mu,                   Defined after [Eq. 32 Suh]                                  (3 x (M_2+1)) matrix
%
% OUTPUT:
%   * acceleration_mode,    Possible values:    * 1 (i.e., 'Mode 1' aka 'No external acceleration Mode')
%                                               * 2 (i.e., 'Mode 2' aka 'External acceleration Mode')
%
% Notation: k = 1 is the most recent sample, whereas k = 3 is the oldest one
%
% Author: Livio Bisogni
% © 2021 Livio Bisogni
%_______________________________________________________________________________________________________

    global M_2 gamma;
%     global estimateExtAccCov_COUNTER    % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!

    % Check number of arguments
    narginchk(2,2);
    
    if (~isequal(size(lambda), [3 (M_2 + 1)]))
        error('lambda must be a (3 x (M_2 + 1)) matrix.');
    end
    if (~isequal(size(mu), [3 (M_2 + 1)]))
        error('mu must be a (3 x (M_2 + 1)) matrix.');
    end
        
    
    % [Eq. pre-34 Suh]
    for j = (M_2 + 1):-1:1
        max_temp = -inf;
        for i = 1:3
            max_temp = max(lambda(i,j) - mu(i,j), max_temp);
        end
        if (max_temp >= gamma)
%             estimateExtAccCov_COUNTER = estimateExtAccCov_COUNTER + 1;    % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
            mode = 2;                   % External acceleration detected
%             max_temp    % DEBUG
%             gamma     % DEBUG
            return;
        end
    end
    
    mode = 1;                           % No external acceleration detected

end