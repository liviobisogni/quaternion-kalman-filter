function alpha_wrapped = wrapTo90(alpha)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [https://it.mathworks.com/matlabcentral/answers/324032-how-to-wrap-angle-in-radians-to-pi-2-pi-2]
% WRAPTO90 Wraps given angle to [-90, +90].
%
% INPUT:
%   * alpha,            Angle                                                   scalar      [deg]
%
% OUTPUT:
%   * alpha_wrapped,    Wrapped angle                                           scalar      [deg]
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);
    
%     if (~isscalar(alpha))
%         error('alpha must be a scalar.');
%     end

    % https://it.mathworks.com/matlabcentral/answers/324032-how-to-wrap-angle-in-radians-to-pi-2-pi-2
    tmp = mod(alpha + 90, 180);
    alpha_wrapped = tmp + 180 * (alpha > 0 & tmp == 0) - 90;
    
    % ALTERNATIVE (for range: [-90, 90) )
%     alpha_wrapped = mod(lambda+90,180)-90;       % For range: [-90, 90)
 
end
