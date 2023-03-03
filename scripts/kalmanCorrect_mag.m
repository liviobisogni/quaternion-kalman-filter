%% ************************************************************************
% ***********    MAGNETOMETER MEASUREMENT UPDATE ("CORRECT")    ***********
% ***********               Author: Livio Bisogni               ***********
% *************************************************************************
%
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               INSTRUCTIONS
%__________________________________________________________________________
%
% Please refer to:  * Suh - Section III. TWO-STEP MEASUREMENT UPDATES
%                   * Trawny - Section 3.2 (Kalman Filter Update).
% It adjusts the projected estimate by an actual magnetometer measurement
% at that time.
%
% Notation used:
%           prev = k-1  (or k, respectively)
%           next = k    (or k+1, respectively)
%           e.g.: q_next_prev denotes q(k|k-1) (or q(k+1|k), respectively)
%__________________________________________________________________________


% function [q_hat__next_next, x_hat__next_next, P__next_next] = kalmanCorrect_mag(...
function [x_hat__next_next, P__next_next] = kalmanCorrect_mag(...
    y_m__next, q_hat__next_next, x_hat_a__next_next, P_a__next_next)

global R_m
global m_tilde


%% 0. Compute estimated rotation matrix
% [Eq. 1 Suh]
% C_hat_n_b: Estimated rotation matrix
C_hat_n_b = quat2RotMat(q_hat__next_next);



%% 1. Compute magnetometer measurement matrix:
% [Eq. 22b Suh]
% H_m: Magnetometer measurement matrix
H_m = [2 * skewSymmetric(C_hat_n_b * m_tilde), zeros(3), zeros(3)];



%% 2. Compute magnetometer measurement error according to:
% [Eq. 22c Suh]
% Magnetometer measurement error
z_m = y_m__next - C_hat_n_b * m_tilde;



%% 3. Update the covariance matrix:
% Suh proposes to apply this update only to the quaternion and not to the
% gyroscope bias by simplifying the state covariance matrix:
% [Eq. 21a Suh]
P_m__next_prev = [P_a__next_next(1:3,1:3),  zeros(3,6);
                  zeros(6,3),               zeros(6,6)];



%% 4. Compute magnetometer residual covariance matrix:
% [Eq. 21b (partial) Suh]
% S_m: Magnetometer residual covariance matrix
S_m = H_m * P_m__next_prev * H_m' + R_m;



%% 5. Compute the magnetometer Kalman gain:
% To limit the effect of the correction only to the yaw component, Suh
% proposed to modify the computation of the gain as following:
% [Eq. 22a Suh]
r_3 = C_hat_n_b * [0; 0; 1];
% and then:
% [Eq. 21b (partial) Suh]
% Suh_matrix: A fictitious matrix proposed by Suh
Suh_matrix = [r_3 * r_3' 	zeros(3,6);
              zeros(6,3)    zeros(6,6)];

% [Eq. 21b Suh]
% K_m: Magnetometer Kalman gain
K_m = Suh_matrix * P_m__next_prev * H_m' / S_m;



%% 6. Compute residual in the magnetometer measurement update:
% [Eq. 21c (partial) Suh]
% r_m: Residual in the magnetometer measurement update
r_m = z_m - H_m * x_hat_a__next_next;



%% 7. Update the state vector:
% [Eq. 21c Suh]
% x_hat__next_next: Updated state
x_hat__next_next = x_hat_a__next_next + K_m * (r_m);



%% 8. Update the covariance matrix:
% [Eq. 21d Suh]
% P__next_next: State covariance matrix
temp = eye(9) - K_m * H_m;
P__next_next = temp * P_a__next_next * temp' + K_m * R_m * K_m';


end