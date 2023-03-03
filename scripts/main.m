% *************************************************************************
% ***********    ATTITUDE ESTIMATION USING A QUATERNION-BASED   ***********
% ***********     KALMAN FILTER WITH ADAPTIVE AND NORM-BASED    ***********
% ***********        ESTIMATION OF EXTERNAL ACCELERATION        ***********
% ***********                                                   ***********
% ***********               Author: Livio Bisogni               ***********
% *************************************************************************
%
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               INSTRUCTIONS
%__________________________________________________________________________
%
% Please refer to:
%            1) Y. S. Suh, "Orientation Estimation Using a Quaternion-Based
%               Indirect Kalman Filter With Adaptive Estimation of External
%               Acceleration," in IEEE Transactions on Instrumentation and
%               Measurement, vol. 59, no. 12, pp. 3296-3305, Dec. 2010,
%               doi: 10.1109/TIM.2010.2047157.
%            2) Nez, A.; Fradet, L.; Marin, F.; Monnet, T.; Lacouture, P.
%               "Identification of Noise Covariance Matrices to Improve
%               Orientation Estimation by Kalman Filter." Sensors 2018, 18,
%               3490. https://doi.org/10.3390/s18103490
%            3) Trawny, N.; Roumeliotis, S.I., "Indirect Kalman filter for
%               3D attitude estimation." In Technical Report; Department of
%               Computer Science and Engineering, University of Minnesota:
%               Minneapolis, MN, USA, 2005.
%            4) Sabatini, A. M. Quaternion-based extended Kalman filter for
%               determining orientation by inertial and magnetic sensing.
%               IEEE Trans. Biomed. Eng., vol. 53, no. 7, 1346–1356, 
%               Jul. 2006.
%            5) Schneider, R; Georgakis, C. How To NOT Make the Extended
%               Kalman Filter Fail. Industrial & Engineering Chemistry
%               Research 2013 52 (9), 3354-3362. doi:10.1021/ie300415d
%            6) Bisogni, L.; Mollaiyan, R.; Pettinari, M.; Neri, P.;
%               Gabiccini, M. Automatic Calibration of a Two-Axis Rotary
%               Table for 3D Scanning Purposes. Sensors 2020, 20, 7107.
%               doi:10.3390/s20247107
%
% Notation used:
%           prev = k-1  (or k, respectively)
%           next = k    (or k+1, respectively)
%           e.g.: q_next_prev denotes q(k|k-1) (or q(k+1|k), respectively)
%
% Please select one of the following sensor_data types:
%       * 'synthetic':  Computer-generated sequence of angular velocity,
%                       acceleration, and magnetic field; see
%                       'generate_synthetic_data.m' for further details.
%       * 'real':       Load sequence of angular velocity, acceleration,
%                       magnetic field, and attitude, measured by real
%                       sensors (gyroscope, accelerometer, magnetometer,
%                       and inclinometer, respectively); see
%                       'import_real_data.m' for further details.
%
%
% close all; clear; clc
sensor_data = 'synthetic';
%__________________________________________________________________________

global N_samples N_samples__theoretical



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                         DATA GENERATION / LOAD
%%_________________________________________________________________________

if strcmpi(sensor_data, 'synthetic')
    % Computer-generated true (reference) and measurement data
    generate_synthetic_data
    N_samp = N_samples;
elseif strcmpi(sensor_data, 'real')
    % Load measurement data from real sensors
    import_real_data
    if (interpolation == 1 || interpolation == 2 || interpolation == 3)
        N_samp = N_samples__theoretical;
    else
        N_samp = N_samples;
    end
else
    warning('Please select a valid data source. Wrong value for sensor_data. Possible values: ''synthetic'' or ''real''.')
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                 Preallocate memory for variables to be saved
%%_________________________________________________________________________
%
% PREDICTED (estimated)
x_hat__next_prev_TEMP     = zeros(9, N_samp);    % state (est.)
b_hat_g__next_prev_TEMP   = zeros(3, N_samp);    % giro bias (est.)
b_hat_a__next_prev_TEMP   = zeros(3, N_samp);    % accelerometer bias (est.)
q_hat__next_prev_TEMP     = zeros(4, N_samp);    % quaternion (est.) 
P__next_prev_TEMP         = zeros(9, 9, N_samp); % P (est.)
C_hat_n_b__next_prev_TEMP = zeros(3, 3, N_samp); % C_n_b (est.)
w_hat__next_prev_TEMP     = zeros(3, N_samp);    % w (est.)
a_hat__next_prev_TEMP     = zeros(3, N_samp);    % a (est.)
m_hat__next_prev_TEMP     = zeros(3, N_samp);    % m (est.)
%
% CORRECTED (updated)
x_hat__next_next_TEMP     = zeros(9, N_samp);    % state (updated)
b_hat_g__next_next_TEMP   = zeros(3, N_samp);    % giro bias (updated)
b_hat_a__next_next_TEMP   = zeros(3, N_samp);    % accelerometer bias (updated)
q_hat__next_next_TEMP     = zeros(4, N_samp);    % quaternion (updated)
P__next_next_TEMP         = zeros(9, 9, N_samp); % P (updated)
C_hat_n_b__next_next_TEMP = zeros(3, 3, N_samp); % C_n_b (updated)
w_hat__next_next_TEMP     = zeros(3, N_samp);    % w (updated)
a_hat__next_next_TEMP     = zeros(3, N_samp);    % a (updated)
m_hat__next_next_TEMP     = zeros(3, N_samp);    % m (updated)
%
r_a = zeros(3,M_1);
lambda = zeros(3,M_1);
mu = zeros(3,M_1);
%
Q_hat_a_b_SAB_TEMP = zeros(3, 3, N_samp);        % estimated Q_a_b by Sabatini
Q_hat_a_b_SUH_TEMP = zeros(3, 3, N_samp);        % estimated Q_a_b by Suh
%
P__next_next_TEMP(:,:,1) = P__prev_prev;
y_g__prev = y_g(:,1);



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                          ITERATE N_samp TIMES
%%_________________________________________________________________________
for i = 2:N_samp


    %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    %                    Acquire current (k=i) sensor output
    %%_________________________________________________________________________
    y_g__next = y_g(:,i);       % current gyroscope measurement                     [rad / s]
    y_a__next = y_a(:,i);       % current accelerometer measurement                 [m / s^2]
    y_m__next = y_m(:,i);       % current magnetometer measurement                  [microT]


    %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    %                       PREDICTION (using GYROSCOPE)
    %%_________________________________________________________________________
    [q_hat__next_prev, x_hat__next_prev, P__next_prev] = kalmanPredict(...
        y_g__prev, y_g__next, q_hat__prev_prev, x_hat__prev_prev, P__prev_prev);


    %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    %                       1st CORRECTION (using ACCELEROMETER)
    %%_________________________________________________________________________
    [q_hat__next_next, x_hat_a__next_next, P_a__next_next, r_a, lambda, mu, Q_hat_a_b_SAB, Q_hat_a_b_SUH] = kalmanCorrect_acc(...
        y_a__next, q_hat__next_prev, x_hat__next_prev, P__next_prev, r_a, lambda, mu, i);

    
    %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    %                       2nd CORRECTION (using MAGNETOMETER)
    %%_________________________________________________________________________
%     [q_hat__next_next, x_hat__next_next, P__next_next] = kalmanCorrect_mag(...
    [x_hat__next_next, P__next_next] = kalmanCorrect_mag(...
        y_m__next, q_hat__next_next, x_hat_a__next_next, P_a__next_next);

    
    %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    %                       Save current (k = i) variables
    %%_________________________________________________________________________
    %
    % PREDICTED (estimated)
    x_hat__next_prev_TEMP(:,i) = x_hat__next_prev;                              % state (est.)
    b_hat_g__next_prev_TEMP(:,i) = x_hat__next_prev_TEMP(4:6,i);                % giro bias (est.)
    b_hat_a__next_prev_TEMP(:,i) = x_hat__next_prev_TEMP(7:9,i);                % accelerometer bias (est.)
    q_hat__next_prev_TEMP(:,i) = q_hat__next_prev;                              % quaternion (est.)    
    P__next_prev_TEMP(:,:,i) = P__next_prev;                                    % P (est.)
    C_hat_n_b__next_prev_TEMP(:,:,i) = quat2RotMat(q_hat__next_prev_TEMP(:,i)); % C_n_b (est.)
    w_hat__next_prev_TEMP(:,i) = y_g__next - b_hat_g__next_prev_TEMP(:,i);      % w (est.)
%     a_hat__next_prev_TEMP(:,i) = C_hat_n_b__next_prev_TEMP(:,:,i) * g_tilde;    % a (est.)
    a_hat__next_prev_TEMP(:,i) = C_hat_n_b__next_prev_TEMP(:,:,i) * g_tilde - b_hat_a__next_prev_TEMP(:,i); % a (est.)
    m_hat__next_prev_TEMP(:,i) = C_hat_n_b__next_prev_TEMP(:,:,i) * m_tilde;    % m (est.)
    %
    % CORRECTED (updated)
    x_hat__next_next_TEMP(:,i) = x_hat__next_next;                              % state (updated)
    b_hat_g__next_next_TEMP(:,i) = x_hat__next_next_TEMP(4:6,i);                % giro bias (updated)
    b_hat_a__next_next_TEMP(:,i) = x_hat__next_next_TEMP(7:9,i);                % accelerometer bias (updated)
    q_hat__next_next_TEMP(:,i) = q_hat__next_next;                              % quaternion (updated)
    P__next_next_TEMP(:,:,i) = P__next_next;                                    % P (updated)
    C_hat_n_b__next_next_TEMP(:,:,i) = quat2RotMat(q_hat__next_next_TEMP(:,i)); % C_n_b (updated)
    w_hat__next_next_TEMP(:,i) = y_g__next - b_hat_g__next_next_TEMP(:,i);      % w (updated)
%     a_hat__next_next_TEMP(:,i) = C_hat_n_b__next_next_TEMP(:,:,i) * g_tilde;    % a (updated)
    a_hat__next_next_TEMP(:,i) = C_hat_n_b__next_next_TEMP(:,:,i) * g_tilde - b_hat_a__next_next_TEMP(:,i); % a (updated)
%         a_hat__next_next_TEMP(:,i) = y_a - C_hat_n_b__next_next_TEMP(:,:,i) * g_tilde - b_hat_a__next_next_TEMP(:,i); % a (updated)
    m_hat__next_next_TEMP(:,i) = C_hat_n_b__next_next_TEMP(:,:,i) * m_tilde;    % m (updated)
    %
    Q_hat_a_b_SAB_TEMP(:,:,i) = Q_hat_a_b_SAB;
    Q_hat_a_b_SUH_TEMP(:,:,i) = Q_hat_a_b_SUH;

    
    %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
    %                      Increment time istant (k = k + 1)
    %%_________________________________________________________________________
    % Increment time instant: k++
%     time_old = time_new;
    q_hat__prev_prev = q_hat__next_next;
    x_hat__prev_prev = x_hat__next_next;
    P__prev_prev = P__next_next;
    y_g__prev = y_g__next;
    
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                       Save all variables (N_samp)
%%_________________________________________________________________________
%
% PREDICTED (estimated)
x_hat__next_prev = x_hat__next_prev_TEMP;
b_hat_g__next_prev = b_hat_g__next_prev_TEMP;
b_hat_a__next_prev = b_hat_a__next_prev_TEMP;
q_hat__next_prev = q_hat__next_prev_TEMP;
P__next_prev = P__next_prev_TEMP;
C_hat_n_b__next_prev = C_hat_n_b__next_prev_TEMP;
w_hat__next_prev = w_hat__next_prev_TEMP;
a_hat__next_prev = a_hat__next_prev_TEMP;
m_hat__next_prev = m_hat__next_prev_TEMP;
%
% CORRECTED (updated)
x_hat__next_next = x_hat__next_next_TEMP;
b_hat_g__next_next = b_hat_g__next_next_TEMP;
b_hat_a__next_next = b_hat_a__next_next_TEMP;
q_hat__next_next = q_hat__next_next_TEMP;
P__next_next = P__next_next_TEMP;
C_hat_n_b__next_next = C_hat_n_b__next_next_TEMP;
w_hat__next_next = w_hat__next_next_TEMP;
a_hat__next_next = a_hat__next_next_TEMP;
m_hat__next_next = m_hat__next_next_TEMP;
%
% Q Matrices
Q_hat_a_b_SAB = Q_hat_a_b_SAB_TEMP;
Q_hat_a_b_SUH = Q_hat_a_b_SUH_TEMP;
%
clear x_hat__next_prev_TEMP, clear b_hat_g__next_prev_TEMP, clear b_hat_a__next_prev_TEMP, clear q_hat__next_prev_TEMP, clear P__next_prev_TEMP, clear C_hat_n_b__next_prev_TEMP, clear w_hat__next_prev_TEMP, clear a_hat__next_prev_TEMP, clear m_hat__next_prev_TEMP
clear x_hat__next_next_TEMP, clear b_hat_g__next_next_TEMP, clear b_hat_a__next_next_TEMP, clear q_hat__next_next_TEMP, clear P__next_next_TEMP, clear C_hat_n_b__next_next_TEMP, clear w_hat__next_next_TEMP, clear a_hat__next_next_TEMP, clear m_hat__next_next_TEMP
clear Q_hat_a_b_SAB_TEMP, clear Q_hat_a_b_SUH_TEMP



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                           Compute Euler Angles
%%_________________________________________________________________________

%% True Euler Angles
if strcmpi(sensor_data, 'synthetic')
    roll_true  = zeros(1,N_samp); % x
    pitch_true = zeros(1,N_samp); % y
    yaw_true   = zeros(1,N_samp); % z
    
    for i=1:N_samp
        %     r = quat2euler(q(:,i));
        r = RotMat2euler(C_n_b(:,:,i));
        %     r = RotMat2euler(quat2RotMat(q(:,i)));
        roll_true(i)   = r(1);
        pitch_true(i)  = r(2);
        yaw_true(i)    = r(3);
        
        roll_true(i)   = wrapToPi(roll_true(i));                % TEMP (for wrapping Euler Angles)
        pitch_true(i)  = wrapToPi2(pitch_true(i));              % TEMP (for wrapping Euler Angles)
        yaw_true(i)    = wrapToPi(yaw_true(i));                 % TEMP (for wrapping Euler Angles)
    end
elseif strcmpi(sensor_data, 'real')
	% Euler angles already computed by vehicle on-board Kalman Filter: they
	% are taken as Reference
else
    warning('Wrong value for sensor_data. Possible values: ''synthetic'' or ''real''.')
end


%% "Measured" Euler Angles (estimated-from-measurement)
%
roll_meas  = zeros(1,N_samp); % x
pitch_meas = zeros(1,N_samp); % y
yaw_meas   = zeros(1,N_samp); % z
%
for i=1:N_samp
%     r = quat2euler(q_m(:,i));
%     r = RotMat2euler(C_n_b_m(:,:,i));
    r = RotMat2euler(quat2RotMat(q_m(:,i)));
    if strcmpi(sensor_data, 'synthetic')
        roll_meas(i)   = r(1);
        pitch_meas(i)  = r(2);
        yaw_meas(i)    = r(3);
    elseif strcmpi(sensor_data, 'real')
        roll_meas(i)   = r(1) + pi;
        pitch_meas(i)  = -r(2);
        yaw_meas(i)    = -r(3);
%
%         roll_meas(i)   = -r(1); % TEMP
%         pitch_meas(i)  = -r(2); % TEMP
%         yaw_meas(i)    = r(3);  % TEMP
    end
    
    roll_meas(i)   = wrapToPi(roll_meas(i));                    % TEMP (for wrapping Euler Angles)
    pitch_meas(i)  = wrapToPi2(pitch_meas(i));                  % TEMP (for wrapping Euler Angles)
    yaw_meas(i)    = wrapToPi(yaw_meas(i));                     % TEMP (for wrapping Euler Angles)
end


%% Estimated Euler Angles
%
roll_hat  = zeros(1,N_samp); % x
pitch_hat = zeros(1,N_samp); % y
yaw_hat   = zeros(1,N_samp); % z
%
for i=1:N_samp
%     r = quat2euler(q_hat(:,i));
    r = RotMat2euler(C_hat_n_b__next_next(:,:,i));
%     r = RotMat2euler(quat2RotMat(q_hat__next_next(:,i)));
    if strcmpi(sensor_data, 'synthetic')
        roll_hat(i)   = r(1);
        pitch_hat(i)  = r(2);
        yaw_hat(i)    = r(3);
    elseif strcmpi(sensor_data, 'real')
        roll_hat(i)   = r(1) + pi;
        pitch_hat(i)  = -r(2);
        yaw_hat(i)    = -r(3);
        
%         roll_hat(i)   = -r(1);  % TEMP
%         pitch_hat(i)  = -r(2);  % TEMP
%         yaw_hat(i)    = r(3);   % TEMP
    end
    
    roll_hat(i)   = wrapToPi(roll_hat(i));                      % TEMP (for wrapping Euler Angles)
    pitch_hat(i)  = wrapToPi2(pitch_hat(i));                    % TEMP (for wrapping Euler Angles)
    yaw_hat(i)    = wrapToPi(yaw_hat(i));                       % TEMP (for wrapping Euler Angles)
end


%% Euler Angles, from Accelerometer & Magnetometer measurements
% [https://habr.com/en/post/499190/]
%
% Estimated Roll (from Accelerometer)
roll_accmag  =  (atan2(a_hat__next_next(2,:), a_hat__next_next(3,:)));                                              % [rad]
%
% Estimated Pitch (from Accelerometer)
pitch_accmag = -wrapToPi(atan2(a_hat__next_next(1,:), sqrt(a_hat__next_next(2,:).^2 + a_hat__next_next(3,:).^2)));  % [rad]
%
% Estimated Yaw (from Magnetometer)
k_0 = 2;
C_n_b_mag = zeros(3,3,N_samp);
C_b_n_mag = zeros(3,3,N_samp);
m_hat__next_next_NEW = zeros(3,N_samp);
yaw_accmag = zeros(1,N_samp);                                                                                      % [rad]
for i=k_0:N_samp
    C_n_b_mag(:,:,i) = euler2RotMat([roll_accmag(i); pitch_accmag(i); 0]);
%     C_n_b_mag(:,:,i) = rotZ(0) * rotY(pitch_accmag(i)) * rotX(roll_accmag(i));                      % CANC !!!!!!!!!!!!!!!
    C_b_n_mag(:,:,i) = C_n_b_mag(:,:,i)';
    m_hat__next_next_NEW(:,i) = C_b_n_mag(:,:,i) * m_hat__next_next(:,i);
    yaw_accmag(i)   =  - atan2(m_hat__next_next_NEW(2,i), m_hat__next_next_NEW(1,i));
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                RMS Errors
%%_________________________________________________________________________
% if strcmpi(sensor_data, 'synthetic')
    
    %%% Estimated - True
    roll_error_hat = roll_hat - roll_true;
    pitch_error_hat = pitch_hat - pitch_true;
    yaw_error_hat = yaw_hat - yaw_true;
    %
    roll_error_hat = wrapTo180(rad2deg(roll_error_hat));
    pitch_error_hat = wrapTo90(rad2deg(pitch_error_hat));
    yaw_error_hat = wrapTo180(rad2deg(yaw_error_hat));
    rms_roll_hat = sqrt(sum(roll_error_hat.^2, 2) / (N_samp - 1));
    rms_pitch_hat = sqrt(sum(pitch_error_hat.^2, 2) / (N_samp - 1));
    rms_yaw_hat = sqrt(sum(yaw_error_hat.^2, 2) / (N_samp - 1));
    disp(['RMSE_roll = ', num2str(rms_roll_hat),  '°,', 9, 'RMSE_pitch = ', num2str(rms_pitch_hat), '°,', 9, 'RMSE_yaw = ', num2str(rms_yaw_hat), '°', 9, '(Kalman)'])

    %%% Measured - True
    roll_error_meas = roll_meas - roll_true;
    pitch_error_meas = pitch_meas - pitch_true;
    yaw_error_meas = yaw_meas - yaw_true;
    %
    roll_error_meas = wrapTo180(rad2deg(roll_error_meas));
    pitch_error_meas = wrapTo90(rad2deg(pitch_error_meas));
    yaw_error_meas = wrapTo180(rad2deg(yaw_error_meas));
    rms_roll_meas = sqrt(sum(roll_error_meas.^2, 2) / (N_samp - 1));
    rms_pitch_meas = sqrt(sum(pitch_error_meas.^2, 2) / (N_samp - 1));
    rms_yaw_meas = sqrt(sum(yaw_error_meas.^2, 2) / (N_samp - 1));
    disp(['RMSE_roll = ', num2str(rms_roll_meas),  '°,', 9, 'RMSE_pitch = ', num2str(rms_pitch_meas), '°,', 9, 'RMSE_yaw = ', num2str(rms_yaw_meas), '°', 9, '(Measurements)'])
    
% elseif strcmpi(sensor_data, 'real')
%     warning('No RMS errors computed, since real sensor data have been used.')
% else
%     warning('Wrong value for sensor_data. Possible values: ''synthetic'' or ''real''.')
% end



% %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% %                            RMS Errors NEW (Sabatini)
% %                                (does NOT work)
% %%_________________________________________________________________________
% if strcmpi(sensor_data, 'synthetic')
%     
%     %%% Estimated - True
%     delta_q_hat = zeros(4,N_samp);
%     delta_theta_hat = zeros(1,N_samp);
%     for i=1:N_samp
%         delta_q_hat(:,i) = quatMultiplication(quatInverse(q(:,i)),q_hat__next_next(:,i));
% %         delta_q_hat(:,i) = q_hat__next_next(:,i);                                           % CANC
%         delta_theta_hat(i) = 2 * acos(delta_q_hat(1,i));
%         delta_theta_hat(i) = wrapTo180(rad2deg(delta_theta_hat(i)));
%     end
%     rms_theta_hat = sqrt(sum(delta_theta_hat.^2, 2) / (N_samp - 1));         % new
%     disp(['RMSE_theta = ', num2str(rms_theta_hat), '°', 9, '(Kalman)'])
%     
%     
%     %%% Measured - True
%     delta_q_meas = zeros(4,N_samp);
%     delta_theta_meas = zeros(1,N_samp);
%     for i=1:N_samp
%         delta_q_meas(:,i) = quatMultiplication(quatInverse(q(:,i)),q_m(:,i));
% %         delta_q_meas(:,i) = q_m(:,i);                                           % CANC
%         delta_theta_meas(i) = 2 * acos(delta_q_meas(1,i));
%         delta_theta_meas(i) = wrapTo180(rad2deg(delta_theta_meas(i)));
%     end
%     rms_theta_meas = sqrt(sum(delta_theta_meas.^2, 2) / (N_samp - 1));         % new
%     disp(['RMSE_theta = ', num2str(rms_theta_meas), '°', 9, '(Measurements)'])
% 
% elseif strcmpi(sensor_data, 'real')
%     warning('No RMS errors computed, since real sensor data have been used.')
% else
%     warning('Wrong value for sensor_data. Possible values: ''synthetic'' or ''real''.')
% end



% %% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% %                               BIASES ESTIMATION:
% %                              MEANS AND VARIANCES
% %%_________________________________________________________________________
% if strcmpi(sensor_data, 'synthetic')
% 
% %     format longEng
% 
%     % Gyro Bias
%     disp(' ')
%     b_g1_error = b_hat_g__next_next(1,:) - b_g(1,:);
%     b_g1_variance = sum(b_g1_error.^2, 2) / (N_samp - 1);
%     b_g1_average = mean(b_g1_error);
%     disp(['b_hat_g,1 - b_g,1:', 9,'VARIANCE(b_g1) = ', num2str(b_g1_variance),  ' [rad^2 / s^2],', 9, 'AVERAGE(b_g1) = ', num2str(b_g1_average), ' [rad / s]'])
%     %
%     b_g2_error = b_hat_g__next_next(2,:) - b_g(2,:);
%     b_g2_variance = sum(b_g2_error.^2, 2) / (N_samp - 1);
%     b_g2_average = mean(b_g2_error);
%     disp(['b_hat_g,2 - b_g,2:', 9,'VARIANCE(b_g2) = ', num2str(b_g2_variance),  ' [rad^2 / s^2],', 9, 'AVERAGE(b_g2) = ', num2str(b_g2_average), ' [rad / s]'])
%     %
%     b_g3_error = b_hat_g__next_next(3,:) - b_g(3,:);
%     b_g3_variance = sum(b_g3_error.^2, 2) / (N_samp - 1);
%     b_g3_average = mean(b_g3_error);
%     disp(['b_hat_g,3 - b_g,3:', 9,'VARIANCE(b_g3) = ', num2str(b_g3_variance),  ' [rad^2 / s^2],', 9, 'AVERAGE(b_g3) = ', num2str(b_g3_average), ' [rad / s]'])
% 
%     % Accelerometer Bias
%     disp(' ')
%     b_a1_error = b_hat_a__next_next(1,:) - b_a(1,:);
%     b_a1_variance = sum(b_a1_error.^2, 2) / (N_samp - 1);
%     b_a1_average = mean(b_a1_error);
%     disp(['b_hat_a,1 - b_a,1:', 9,'VARIANCE(b_a1) = ', num2str(b_a1_variance),  ' [rad^2 / s^2],', 9, 'AVERAGE(b_a1) = ', num2str(b_a1_average), ' [rad / s]'])
%     %
%     b_a2_error = b_hat_a__next_next(2,:) - b_a(2,:);
%     b_a2_variance = sum(b_a2_error.^2, 2) / (N_samp - 1);
%     b_a2_average = mean(b_a2_error);
%     disp(['b_hat_a,2 - b_a,2:', 9,'VARIANCE(b_a2) = ', num2str(b_a2_variance),  ' [rad^2 / s^2],', 9, 'AVERAGE(b_a2) = ', num2str(b_a2_average), ' [rad / s]'])
%     %
%     b_a3_error = b_hat_a__next_next(3,:) - b_a(3,:);
%     b_a3_variance = sum(b_a3_error.^2, 2) / (N_samp - 1);
%     b_a3_average = mean(b_a3_error);
%     disp(['b_hat_a,3 - b_a,3:', 9,'VARIANCE(b_a3) = ', num2str(b_a3_variance),  ' [rad^2 / s^2],', 9, 'AVERAGE(b_a3) = ', num2str(b_a3_average), ' [rad / s]'])
% 
% elseif strcmpi(sensor_data, 'real')
%     warning('No RMS errors computed, since real sensor data have been used.')
% else
%     warning('Wrong value for sensor_data. Possible values: ''synthetic'' or ''real''.')
% end



%% Counter of 'External acceleration detected' (per function):
% N_samp                                  % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
% estimateExtAccCov_COUNTER               % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
% acc_norm_adapt_COUNTER                  % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                              Save workspace
%%_________________________________________________________________________
% cd('/Users/v/Documents/MATLAB/attitude_estimation/images')
% save('KF_bisogni_MODIFYDATE.mat')
% cd('/Users/v/Documents/MATLAB/attitude_estimation')



% counter_instants_above_threshold                                % DEBUGGING PRINT
%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               Plot results
%%_________________________________________________________________________
plot_results


% save('workspace_synth_200Hz_Suh_multconst_20210606.mat');