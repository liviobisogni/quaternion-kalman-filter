% *************************************************************************
% ***********              TIME UPDATE ("PREDICT")              ***********
% ***********               Author: Livio Bisogni               ***********
% *************************************************************************
%
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               INSTRUCTIONS
%__________________________________________________________________________
%
% Please refer to:  * Suh - Section II. QUATERNION-BASED INDIRECT KALMAN FILTER
%                   * Trawny - Section 2.6 (Propagation Equations).
% It projects the current state estimate ahead in time by two actual
% gyroscope measurements: at that time and the previous one.
%
% Notation used:
%           prev = k-1  (or k, respectively)
%           next = k    (or k+1, respectively)
%           e.g.: q_next_prev denotes q(k|k-1) (or q(k+1|k), respectively)
%__________________________________________________________________________


function [q_hat__next_prev, x_hat__next_prev, P__next_prev] = kalmanPredict(...
    y_g__prev, y_g__next, q_hat__prev_prev, x_hat__prev_prev, P__prev_prev)


%% 1. Propagate the quaternion using a third-order integrator

% [Eq. 18 Suh]
q_hat__next_prev = quatFirstIntegration(q_hat__prev_prev, y_g__next, y_g__prev, 'Suh');   % Possible values: 'Suh' or 'Trawny' or 'Yuan'



%% 2. Compute the state transition matrix Φ and the discrete-time process noise covariance matrix Qd.

% [Eq. 16a Suh], 'approximated'; OR [Eq. post-15 Suh], 'precise'; OR [Eq. 187 Trawny], 'Trawny')
% Phi: State transition matrix
Phi = Phi_matrix(y_g__next, 'approximated');    % 'precise' OR 'approximated' OR 'Trawny' (('Trawny' actually DOES NOT work))
% NB: Phi_matrix.m 'precise' è problematico in presenza di accelerazioni esterne && di estimateExtAccCov.m !!!!!!!!!!

% [Eq. 16b Suh], 'approximated'; OR [Eq. post-15 Suh], 'precise'; OR [Eq. 208 Trawny], 'Trawny')
% Q_d: Noise covariance Matrix
Q_d = Q_d_matrix(y_g__next, 'approximated');    % 'precise' OR 'approximated' OR 'Trawny' (('Trawny' actually DOES NOT work))



%% 3. Project the state x and the state covariance matrix P ahead

% [Eq. 17a Suh]
% Project the state ahead
x_hat__next_prev = Phi * x_hat__prev_prev;

% [Eq. 17b Suh]
% Project the state covariance matrix ahead
P__next_prev = Phi * P__prev_prev * Phi' + Q_d;


end