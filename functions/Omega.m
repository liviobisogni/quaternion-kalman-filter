function Omega_omega = Omega(omega)
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% [Eq. A12 ("A Kalman Filter for Nonlinear Attitude Estimation Using Time Variable Matrices and
% Quaternions" - Alvaro Deibe)]
% OMEGA Computes the Omega(omega) matrix.
% NB: Omega(omega) appears in the product of a vector and a quaternion, and is used for example in the
%     quaternion derivative.
% (( NB: [Eq. 64 Trawny] DOES NOT work properly anymore, probably since q(1) = q_0 is the scalar part ))
%
% INPUT:
%	* omega_next,       Angular velocity omega(k)                               (3 x 1) vector      [rad / s]
%
% OUTPUT:
%	* Omega__omega,     Omega matrix                                            (3 x 3) matrix
%
% Author: Livio Bisogni
%_______________________________________________________________________________________________________

    % Check number of arguments
    narginchk(1,1);

    if (~isequal(size(omega), [3 1]))
        error('omega must be a (3 x 1) vector.');
    end

%     % [Eq. 64a Trawny] (DOES NOT work!)
%     Omega_omega = [-skewSymmetric(omega),    omega;
%                    -omega',                  0];
%     
%     w1 = omega(1);
%     w2 = omega(2);
%     w3 = omega(3);
%     
%     % [Eq. 64b Trawny] (DOES NOT work!)
%     Omega_omega = [0,   w3,     -w2,    w1;
%                    -w3, 0,      w1,     w2;
%                    w2,  -w1,    0,      w3;
%                    -w1, -w2,    -w3,    0];
    
    
    % [Eq. A12 ("A Kalman Filter for Nonlinear Attitude Estimation Using Time
    %            Variable Matrices and Quaternions" - Alvaro Deibe)]
    w1 = omega(1);
    w2 = omega(2);
    w3 = omega(3);
    
    Omega_omega = [0,   -w1,    -w2,    -w3;
                   w1,  0,      w3,     -w2;
                   w2,  -w3,    0,      w1;
                   w3,  w2,     -w1,    0];

               
      %% ALTERNATIVE equation (it DOES work, btw  :) )     
%     % [Eq. 13 ("A New Quaternion-Based Kalman Filter for Real-Time Attitude
%                 Estimation Using the Two-Step Geometrically-Intuitive
%                 Correction Algorithm" - Kaiqiang Feng)]
%     Omega_omega = [0,           -omega';
%                    omega,       -skewSymmetric(omega)];
    
end