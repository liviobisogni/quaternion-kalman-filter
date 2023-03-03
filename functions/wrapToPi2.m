function alpha_wrapped = wrapToPi2(alpha)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [https://it.mathworks.com/matlabcentral/answers/324032-how-to-wrap-angle-in-radians-to-pi-2-pi-2]
% WRAPTOPI2 Wraps given angle to [-pi/2, +pi/2].
%
% INPUT:
%   * alpha,            Angle                                                   scalar      [rad]
%
% OUTPUT:
%   * alpha_wrapped,    Wrapped angle                                           scalar      [rad]
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
%     if (~isscalar(alpha))
%         error('alpha must be a scalar.');
%     end

    % https://it.mathworks.com/matlabcentral/answers/324032-how-to-wrap-angle-in-radians-to-pi-2-pi-2
    tmp = mod(alpha + pi/2, pi);
    alpha_wrapped = tmp + pi * (alpha > 0 & tmp == 0) - pi/2;
    
    % ALTERNATIVE (for range: [-pi/2, pi/2) )
%     alpha_wrapped = mod(lambda+pi/2,pi)-pi/2;       % For range: [-pi/2, pi/2)
 
end