%% ************************************************************************
% ***********    ACCELEROMETER MEASUREMENT UPDATE ("CORRECT")   ***********
% ***********               Author: Livio Bisogni               ***********
% *************************************************************************
%
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               INSTRUCTIONS
%__________________________________________________________________________
%
% Please refer to:  * Suh - Section III. TWO-STEP MEASUREMENT UPDATES
%                   * Trawny - Section 3.2 (Kalman Filter Update).
% It adjusts the projected estimate by an actual accelerometer measurement
% at that time.
%
% Notation used:
%           prev = k-1  (or k, respectively)
%           next = k    (or k+1, respectively)
%           e.g.: q_next_prev denotes q(k|k-1) (or q(k+1|k), respectively)
%
% Please select (at line 32) one of the following Q_hat_a_b__algorithm:
%                   * 'suh':        use Suh's algorithm
%                   * 'sabatini':   use Sabatini's algorithm
%                   * 'zeros':      use a 3x3 zeros matrix
%
% e.g.,             Q_hat_a_b__algorithm = 'suh';       <--- go to line 32
%__________________________________________________________________________


function [q_hat__next_next, x_hat_a__next_next, P_a__next_next, r_a, lambda, mu, Q_hat_a_b_SAB, Q_hat_a_b_SUH] = kalmanCorrect_acc(...
    y_a__next, q_hat__next_prev, x_hat__next_prev, P__next_prev, r_a, lambda, mu, i)

Q_hat_a_b__algorithm = 'suh';                           % <------- You can select/change the Q_hat_a_b algorithm here

global R_a
global g_tilde
global ext_acc_detection_time_SUH ext_acc_detection_time_SAB


%% 0. Compute estimated rotation matrix
% [Eq. 1 Suh]
% C_hat_n_b__next_prev: Estimated rotation matrix
C_hat_n_b__next_prev = quat2RotMat(q_hat__next_prev);



%% 1. Compute accelerometer measurement matrix:
% [Eq. 19d Suh]
% H_a: Accelerometer measurement matrix
% (it is the noiseless connection between the state vector and the measurement vector)
H_a = [2 * skewSymmetric(C_hat_n_b__next_prev * g_tilde), zeros(3), eye(3)];



%% 2. Compute accelerometer measurement error according to:
% [Eq. 19e Suh]
% Accelerometer measurement error
z_a = y_a__next - C_hat_n_b__next_prev * g_tilde;



%% 3. Estimate external acceleration covariance matrix:
% Q_hat_a_b: Estimated external acceleration covariance matrix
%
%   * Method proposed by Suh ([Eq. 34 - 35 Suh]):
[Q_hat_a_b_SUH, lambda, mu] = estimateExtAccCov_Suh(r_a, lambda, mu, H_a, P__next_prev, R_a);
%
% Q_hat_a_b_SUH = Q_hat_a_b; % DEBUG ONLY
%   * An ALTERNATIVE method, cited by Suh ([Eq. 37 Suh]):
Q_hat_a_b_SAB = estimateExtAccCov_Sab(y_a__next);       % It DOES WORK as well (ALTERNATIVE method proposed by A. M. Sabatini,
                                                        % “Quaternion-based extended Kalman filter for determining orientation
                                                        % by inertial and magnetic sensing,” IEEE Trans. Biomed. Eng., vol. 53,
                                                        % no. 7, pp. 1346–1356, Jul. 2006.)
%
Q_hat_a_b_ZEROS = zeros(3);
%
% *************************************************************************
if strcmpi(Q_hat_a_b__algorithm, 'suh')
    Q_hat_a_b = Q_hat_a_b_SUH;
elseif  strcmpi(Q_hat_a_b__algorithm, 'sabatini')
    Q_hat_a_b = Q_hat_a_b_SAB;
elseif strcmpi(Q_hat_a_b__algorithm, 'zeros')
    Q_hat_a_b = Q_hat_a_b_ZEROS;
else
    warning('Please select a valid Q_hat_a_b algorithm. Wrong value for Q_hat_a_b__algorithm. Possible values: ''suh'', ''sabatini'' or ''zeros''.')
end
% *************************************************************************
if (~isequal(Q_hat_a_b_SUH ,zeros(3)))
%     i                                       % DEBUG; it prints External acceleration detection instants
    ext_acc_detection_time_SUH(i) = 1;
end
if (~isequal(Q_hat_a_b_SAB ,zeros(3)))
%     i                                       % DEBUG; it prints External acceleration detection instants
    ext_acc_detection_time_SAB(i) = 1;
end


%% 4. Compute accelerometer residual covariance matrix:
% [Eq. 19a (partial) Suh]
% S_a: Accelerometer residual covariance matrix
S_a = H_a * P__next_prev * H_a' + R_a + Q_hat_a_b;



%% 5. Compute accelerometer Kalman gain:
% [Eq. 19a Suh]
% K_a: Accelerometer Kalman gain
K_a = P__next_prev * H_a' / S_a;



% shift each row of r_a 1 position to the right (i.e., age r_a's of 1 sample)
r_a = circshift(r_a,1,2);

%% 6. Compute residual in the accelerometer measurement update:
% [Eq. 19b (partial) Suh and Eq. pre-29 Suh]
% r_a: Residual in the accelerometer measurement update
r_a(:,1) = z_a - H_a * x_hat__next_prev;                      % r_a(:,1) = the most recent r_a



%% 7. Update the state vector:
% [Eq. 19b Suh]
% x_hat_a__next_next: Updated state
x_hat_a__next_next = x_hat__next_prev + K_a * (r_a(:,1));



%% 8. Update the covariance matrix:
% [Eq. 19c Suh]
% P_a__next_next: State covariance matrix
temp = eye(9) - K_a * H_a;
P_a__next_next = temp * P__next_prev * temp' + K_a * (R_a + Q_hat_a_b) * K_a';



%% 9. Update the attitude quaternion

% [Eq. 20a Suh]
% Extract the imaginary part of the attitude error quaternion from the state vector
% q_e : Orientation error (3x1) vector (i.e, the vector part of the error quaternion)
q_e = x_hat_a__next_next(1:3);

% [Eq. 6 Suh]
% Create the attitude error quaternion
% q_tilde_e: Orientation error (4x1) quaternion
q_tilde_e = [1; q_e];

% [Eq. 20b Suh] (and [Eq. 5 Suh])
% Compute the new attitude quaternion
q_hat__next_next = quatMultiplication(q_hat__next_prev, q_tilde_e);

% [Eq. 20c Suh]
% Normalize the quaternion
q_hat__next_next = q_hat__next_next / norm(q_hat__next_next);

% Ensure quaternion scalar part is non-negative:
if (q_hat__next_next(1) < 0)
    q_hat__next_next = -q_hat__next_next;
end

% [Eq. 20d Suh]
% Set q_e to zero (i.e., set the first 3 elements of the state array to 0)
x_hat_a__next_next(1:3) = zeros(3,1);


end