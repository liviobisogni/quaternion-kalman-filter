function x = rand_range(xmin, xmax)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [https://it.mathworks.com/matlabcentral/answers/66763-generate-random-numbers-in-range-from-0-8-to-4]
% RAND_RANGE Generates a value from the uniform distribution on the interval (xmin, xmax).
%
% INPUT:
%   * xmin,             Minimum possible value                                  scalar      []
%   * xmax,             Maximum possible value                                  scalar      []
%
% OUTPUT:
%   * x,                Pseudo-random number in (xmin, xmax)                    scalar      []
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(2,2);
    
    if (~isscalar(xmin))
        error('xmin must be a scalar.');
    end
    if (~isscalar(xmax))
        error('xmax must be a scalar.');
    end
    
    x = xmin + (xmax - xmin) .* rand(1);

end