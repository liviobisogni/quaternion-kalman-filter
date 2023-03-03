%% ************************************************************************
% ***********              IMPORT REAL SENSOR DATA              ***********
% ***********               Author: Livio Bisogni               ***********
% *************************************************************************
%
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               INSTRUCTIONS
%__________________________________________________________________________
%
% It loads the real data from the 'data' folder (which one in particular is
% chosen after in the code).
%
% Notation used:
%           prev = k-1  (or k, respectively)
%           next = k    (or k+1, respectively)
%           e.g.: q_next_prev denotes q(k|k-1) (or q(k+1|k), respectively)
%
% Please select one of the following interpolation modes:
interpolation = 0;  % * if '1', interpolate samples (using Matlab function 'interp1(nearest)')
                    % * if '2', linearly interpolate samples (using 'lin_interpolate', located in the 'functions' directory);
                    %   it is a time-consuming task
                    % * if '3', load previously interpolated data
                    % otherwise do nothing
%
% *** !!! Please Note: paths should be adjusted according to your real data
% folder location. !!! ***
%
%__________________________________________________________________________




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                              GLOBAL VARIABLES
%%_________________________________________________________________________
global dt N_samples N_samples__theoretical         % time-related stuff
global y_g                  % measured gyroscope data
global y_a                  % measured accelerometer data
global y_m                  % measured magnetometer data

global R_a R_m                      % covariance measurement matrices of the accelerometer / magnetometer
global g_standard g_tilde m_tilde
global M_1 M_2 gamma
global Q
global s
global epsilon_a
% global ext_acc_detection_time       % for printing instants where external acceleration is detected by norm-based algorithm
global ext_acc_detection_time_SUH ext_acc_detection_time_SAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG
% Counter of 'External acceleration detected' (per function):
global estimateExtAccCov_COUNTER    % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
global acc_norm_adapt_COUNTER       % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
estimateExtAccCov_COUNTER = 0;      % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
acc_norm_adapt_COUNTER = 0;         % TEMP (DEBUG ONLY)!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                 LOAD
%%_________________________________________________________________________
% Please chose one of them:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (interpolation == 3)
    load('/Users/v/Documents/MATLAB/attitude_estimation/images/Real DATA_interpolated___test_volo_20200213__volo_3__LOG00041_parsed_seg5 Suh/DATA_interpolated___test_volo_20200213__volo_3__LOG00041_parsed_seg5__GLOBAL_WORKSPACE.mat')
    return
end

% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200211/test 1/LOG00035.mat')
% minutes_of_silence = 2;

% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200218/LOG00076_parsed_seg4.mat')
% minutes_of_silence = 2;

load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200213/volo 3/LOG00041_parsed_seg5.mat')   % YEAH
minutes_of_silence = 5.5;           % <--- useful to estimate gyro and accelerometer biases

% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200218/LOG00079_parsed_seg2.mat')
% minutes_of_silence = 0.75;
% 
% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200218/LOG00076_parsed_seg4.mat')
% minutes_of_silence = 6;
% 
% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200225/LOG00102_parsed_seg1.mat')
% minutes_of_silence = 3.5;
% 
% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200225/LOG00101_parsed_seg2.mat') % NO 
% minutes_of_silence = 2.5;
% 
% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200304/LOG00105_parsed_seg2.mat') % NO
% minutes_of_silence = 3.5; % no silence......
% 
% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200320/TestVolo.mat')
% minutes_of_silence = 1;
% 
% load('/Users/v/Documents/MATLAB/attitude_estimation/data/flight_data/test_volo_20200213/volo 1 bad pos/LOG00038_parsed_seg1.mat') % OK, ma accelerazione rumorosa
% minutes_of_silence = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                 CONSTANTS
%%_________________________________________________________________________
% [Eq. post-3 Suh]
% g_standard: Standard acceleration of gravity
g_standard = 9.80665;                                           % [m / s^2]                             <-------------------
% g_standard = norm(mean(DATA(1:minutes_of_silence*60, 8:10))');                                        % !!!!!!!!!!!!!!!!!!!!!!!!!!
% g_tilde: Gravitational field, expressed in the navigation frame
g_tilde = [0; 0; g_standard];

% [Eq. post-3 Suh]
% Data source:      https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
% dip_angle: Dip angle (aka magnetic inclination) @ Viareggio (43° 52' 52" N, 10° 14' 6" E), IT, February 11th, 2021
% P.S.: (43° 52' 52" N, 10° 14' 6" E) = (43.88111° N, 10.23500° E)
dip_angle = deg2rad(60.10803); % 60° 6' 29"                     % [rad]                                 <-------------------
% m_tilde_intensity: Magnetic field intensity
m_tilde_intensity = 47.1179;                                    % [microT]                              <-------------------
% m_tilde_intensity = norm(mean(DATA(1:minutes_of_silence*60, 11:13))') * 100;                          % !!!!!!!!!!!!!!!!!!!!!!!!!!
% m_tilde: Magnetic field (simpler model), expressed in the navigation frame
m_tilde = [cos(dip_angle); 0; -sin(dip_angle)] * m_tilde_intensity;
% m_tilde: Magnetic field (wicked model) predicted @ Viareggio (43° 52' 52" N, 10° 14' 6" E), IT, February 11th, 2021, expressed in the navigation frame
% m_tilde = [23.4485; 1.2545; -40.8497];                          % [microT]                              <-------------------


% Adaptive Suh's algorithm parameters
% [Eq. 36 Suh]
Q_b_g = 0.0000001 * eye(3);                                     % 0.000001                              <------------------- !
Q_b_a = 0.00001 * eye(3);                                       % 0.000001                              <------------------- !
M_1 = 3;                                                        % []                                    <------------------- (not used)
M_2 = 2;                                                        % []                                    <------------------- (not used)
gamma = 0.1;                                                    %                                       <------------------- (not used)

% Q_b_g = 0.0000001 * eye(3);
% Q_b_a = 0.00001 * eye(3);

% Q_b_g = 0.000001 * eye(3);
% Q_b_a = 0.00005 * eye(3);


% Norm-based Sabatini's algorithm parameters (used by Suh)
% [Eq. 37 Suh]
s = 10;                                                         % 10; try also 1 and 100                <-------------------
epsilon_a = 0.25;                                               % [m / s^2]                             <-------------------



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                   TIME
%%_________________________________________________________________________

% N_samples: Number of samples
N_samples = size(T,1);                                          % []

% time: Time stamp array [..., time(k-1), time(k), time(k+1), ...]; it may
% NOT have constant time step!!!
time = T;                                                       % [s]

dt = T(2) - T(1);                                               % [s]
f_sampling = 1 / dt;                                            % [Hz]

% flight_duration: Simulation duration
flight_duration = T(N_samples);                                 % [s]                                   <-------------------

N_samples__theoretical = floor(flight_duration * f_sampling);   % []
if (interpolation == 1 || interpolation == 2)
    N_samp = N_samples__theoretical;
else
    N_samp = N_samples;
end

% ext_acc_detection_time: = 1 only when external acceleration is detected; otherwise (in all other time steps) it remains null
ext_acc_detection_time_SAB = zeros(1,N_samp);
ext_acc_detection_time_SUH = zeros(1,N_samp);

T_theoretical = zeros(N_samp, 1);
for j=1:(N_samp - 0)
    T_theoretical(j) = j * dt;                          % {0.02, 0.04, 0.06, ...}
end

dt_diff_real = zeros(N_samples, 1);
dt_diff = zeros(N_samples, 1);
dt_diff_diff = zeros(N_samples, 1);                     % OLD, CANC??
dt_diff_real(1) = T(1) - 0;
dt_diff(1) = T_theoretical(1) - T(1);
dt_diff_diff(1) = 0;
for j=2:(N_samples - 0)
    dt_diff_real(j) = T(j) - T(j-1);
    dt_diff(j) = abs(T_theoretical(j) - T(j));          % OLD, CANC??
    dt_diff_diff(j) = dt_diff(j) - dt_diff(j-1);        % OLD, CANC??
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                            ESTIMATED EULER ANGLES
%                      *** N.B.: taken as reference! ***
%%_________________________________________________________________________
if (interpolation == 1)
    roll_true  = interp1(T', DATA(:,2)', dt:dt:flight_duration, 'nearest');    % [rad]
    pitch_true = interp1(T', DATA(:,3)', dt:dt:flight_duration, 'nearest');    % [rad]
    yaw_true   = interp1(T', DATA(:,4)', dt:dt:flight_duration, 'nearest');    % [rad]
elseif (interpolation == 2)
    roll_true  = lin_interpolate(T', DATA(:,2)', N_samples, N_samp);    % [rad]
    pitch_true = lin_interpolate(T', DATA(:,3)', N_samples, N_samp);    % [rad]
    yaw_true   = lin_interpolate(T', DATA(:,4)', N_samples, N_samp);    % [rad]
else
    roll_true  = DATA(:,2)';                                        % [rad]
    pitch_true = DATA(:,3)';                                        % [rad]
    yaw_true   = DATA(:,4)';                                        % [rad]
end
% temp = rotX(-pi) * [roll_true; pitch_true; yaw_true];
% roll_true = temp(1,:);
% pitch_true = temp(2,:);
% yaw_true = temp(3,:);


%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                 GYROSCOPE
%%_________________________________________________________________________

% sigma_g: Standard deviation of the gyroscope global errors (noises, external perturbations, calibration defects etc.)
sigma_g = 0.002;  % 0.01;    0.06;   0.001;                     % [rad / s]                             <------------------- !
n_g_variance = sigma_g^2;                                       % [rad^2 / s^2]

% R_g: Covariance measurement matrix of the gyroscope:
R_g = n_g_variance * eye(3);                                    % [rad^2 / s^2]

% y_g: Measured Angular Velocity (i.e., gyro output) [gyro_x; gyro_y; gyro_z]
if (interpolation == 1)
    y_g = zeros(3,N_samp);                                      % [rad / s]
    y_g(1,:) = interp1(T', DATA(:,5)', dt:dt:flight_duration, 'nearest');
    y_g(2,:) = interp1(T', DATA(:,6)', dt:dt:flight_duration, 'nearest');
    y_g(3,:) = interp1(T', DATA(:,7)', dt:dt:flight_duration, 'nearest');
elseif (interpolation == 2)
    y_g = zeros(3,N_samp);                                      % [rad / s]
    y_g(1,:) = lin_interpolate(T', DATA(:,5)', N_samples, N_samp);
    y_g(2,:) = lin_interpolate(T', DATA(:,6)', N_samples, N_samp);
    y_g(3,:) = lin_interpolate(T', DATA(:,7)', N_samples, N_samp);
else
    y_g = DATA(:,5:7)';                                         % [rad / s]
end
% y_g = rotX(-pi) * y_g;
% y_g_1_TEMP = y_g(1,:);
% y_g(1,:) = y_g(2,:);
% y_g(2,:) = y_g_1_TEMP;
% y_g(3,:) = -y_g(3,:);
    
% b_g: Gyroscope bias
b_g = zeros(3, N_samp);
for i=1:N_samp
%     b_g_OLD(:,i) = mean(DATA(1:minutes_of_silence*60, 5:7));
    
    b_g(:,i) = mean(y_g(:, 1:minutes_of_silence*60)')';           % NEW
end




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                 QUATERNION
%%_________________________________________________________________________

% q: true quaternion, i.e., q from w (true omega)
q = zeros(4,N_samp);
% q(:,1) = euler2quat(euler_0);                                                             % one of the first alternatives
% q(:,1) = RotMat2quat(euler2RotMat(euler_0));                                                % used from May 11, 2021 till ...; semi-new!!!
q(:,1) = RotMat2quat(rotZ(yaw_true(1)) * rotX(pitch_true(1)) * rotX(roll_true(1)));                        	% NEW!!!!!!!!!!!!!!!!!!!!!!!!!

% for i=2:N_samp
%     q(:,i) = quatFirstIntegration(q(:,i-1), y_g(:,i), y_g(:,i-1), dt, 'Suh');  % Possible values: 'Suh' or 'Trawny' or 'Yuan'
% end

for i=2:N_samp
%     q(:,i) = RotMat2quat(euler2RotMat([roll_true(i); pitch_true(i); yaw_true(i)]));                                                % used from May 11, 2021 till ...; semi-new!!!
    q(:,i) = RotMat2quat(rotZ(yaw_true(i)) * rotX(pitch_true(i)) * rotX(roll_true(i)));                        	% NEW!!!!!!!!!!!!!!!!!!!!!!!!!
end

% [Eq. 1 Suh] ([Eq. 79 Trawny])
% C_n_b: True associated transformation matrix (function of the true attitude quaternion q)
% C_n_b = zeros(3,3,N_samp);
% for i = 1:N_samp
%     C_n_b(:,:,i) = quat2RotMat(q(:,i));
% end




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               ACCELEROMETER
%%_________________________________________________________________________

% sigma_a: Standard deviation of the accelerometer global errors (noises, external perturbations, calibration defects etc.)
sigma_a = 3; %0.004;     0.039;                                           % [m / s^2]                             <------------------- !
n_a_variance = sigma_a^2;                                       % [m^2 / s^4]

% R_a: Covariance measurement matrix of the accelerometer:
R_a = n_a_variance * eye(3);                                    % [m^2 / s^4]

% y_a  Measured acceleration (i.e., accelerometer output)
if (interpolation == 1)
    y_a = zeros(3,N_samp);                                      % [m / s^2]
    y_a(1,:) = interp1(T', DATA(:,8)', dt:dt:flight_duration, 'nearest');
    y_a(2,:) = interp1(T', DATA(:,9)', dt:dt:flight_duration, 'nearest');
    y_a(3,:) = interp1(T', DATA(:,10)', dt:dt:flight_duration, 'nearest');
elseif (interpolation == 2)
    y_a = zeros(3,N_samp);                                      % [m / s^2]
    y_a(1,:) = lin_interpolate(T', DATA(:,8)', N_samples, N_samp);
    y_a(2,:) = lin_interpolate(T', DATA(:,9)', N_samples, N_samp);
    y_a(3,:) = lin_interpolate(T', DATA(:,10)', N_samples, N_samp);
else
    y_a = DATA(:,8:10)';                                        % [m / s^2]
end
% y_a = rotX(-pi) * y_a;
% y_a_1_TEMP = y_a(1,:);
% y_a(1,:) = y_a(2,:);
% y_a(2,:) = y_a_1_TEMP;
% y_a(3,:) = -y_a(3,:);       % così la gravità da negativa diventa positiva

% b_a: Accelerometer bias
b_a = zeros(3, N_samp);
for i=1:N_samp
%     b_a_OLD(:,i) = mean(DATA(1:minutes_of_silence*60, 8:10))' + g_tilde;
    
    b_a(:,i) = mean(y_a(:, 1:minutes_of_silence*60)')' - g_tilde;           % NEW
    b_a(:,i) = mean(y_a(:, 1:minutes_of_silence*60)')' + g_tilde;           % NEW
end

% % a_b: External acceleration
% a_b = DATA(:,[1 4 7]+15+2)';                % i.e., 18    21    24




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               MAGNETOMETER
%%_________________________________________________________________________

% sigma_m: Standard deviation of the magnetometer global errors (noises, external perturbations, calibration defects etc.)
sigma_m = 8;          %0.2;    [0.022; 0.111; 0.028];           % [microT]                             <------------------- !
n_m_variance = sigma_m^2;                                       % [(microT)^2]

% R_m: Covariance measurement matrix of the magnetometer:
R_m = n_m_variance * eye(3);                                    % [10^-12 T^2]

% y_m: Measured magnetic field (i.e., Magnetometer output)
if (interpolation == 1)
    y_m = zeros(3,N_samp);                                      % [microT]
    y_m(1,:) = interp1(T', DATA(:,11)', dt:dt:flight_duration, 'nearest') * 100;
    y_m(2,:) = interp1(T', DATA(:,12)', dt:dt:flight_duration, 'nearest') * 100;
    y_m(3,:) = interp1(T', DATA(:,13)', dt:dt:flight_duration, 'nearest') * 100;
elseif (interpolation == 2)
    y_m = zeros(3,N_samp);                                      % [microT]
    y_m(1,:) = lin_interpolate(T', DATA(:,11)', N_samples, N_samp) * 100;
    y_m(2,:) = lin_interpolate(T', DATA(:,12)', N_samples, N_samp) * 100;
    y_m(3,:) = -lin_interpolate(T', DATA(:,13)', N_samples, N_samp) * 100;
else
    y_m = DATA(:,11:13)' * 100;                                 % [microT]
end
% y_m = rotX(-pi) * y_m;
% y_m_1_TEMP = y_m(1,:);
% y_m(1,:) = y_m(2,:);
% y_m(2,:) = y_m_1_TEMP;
% y_m(3,:) = -y_m(3,:);       % così la componente z del campo magnetico da positiva diventa negativa




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                     INITIAL ESTIMATED EULER ANGLES AND
%              ESTIMATED QUATERNION BY INTEGRATING ANGULAR VELOCITY
%%_________________________________________________________________________

% Initial Euler Angles, estimated from Accelerometer & Magnetometer measurements:
% [https://habr.com/en/post/499190/]
%
k_0 = 1;                                                        % Index of the first non-null measurement; tipically, either 1 or 2     <-------------------
% Estimated initial Roll (from Accelerometer)
% [Eq. 25 - Tilt Sensing Using a Three-Axis Accelerometer]
roll_accmag_0  = (atan2(y_a(2,k_0), y_a(3,k_0)));                                   % [rad]
%
% Estimated initial Pitch (from Accelerometer)
% [Eq. 26 - Tilt Sensing Using a Three-Axis Accelerometer]
pitch_accmag_0 = -wrapToPi(atan2(y_a(1,k_0), sqrt(y_a(2,k_0).^2 + y_a(3,k_0).^2))); % [rad]
%
% Estimated initial Yaw (from Magnetometer)
% [Eq. 22 - Implementing a Tilt-Compensated eCompass using Accelerometer and Magnetometer Sensors]
% C_n_b__0 = euler2RotMat([roll_accmag_0; pitch_accmag_0; 0]);
% C_n_b__0 = rotZ(0) * rotX(pitch_accmag_0) * rotX(roll_accmag_0);                      % !!!!!!!!!!!!!!! % old, canc?
C_n_b__0 = rotZ(0) * rotY(pitch_accmag_0) * rotX(roll_accmag_0);                      % !!!!!!!!!!!!!!!

C_b_n__0 = C_n_b__0';
y_m_0_NEW = C_b_n__0 * y_m(:,k_0);
yaw_accmag_0 = -atan2(y_m_0_NEW(2), y_m_0_NEW(1));                                  % [rad]

%%%%%%%%%%%%%%% TEMP!!! %%%%%
% roll_accmag_0  = deg2rad(30);
% pitch_accmag_0 = deg2rad(40);
% yaw_accmag_0   = deg2rad(50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimated initial Euler angles
euler_hat_0 = [roll_accmag_0; pitch_accmag_0; yaw_accmag_0];
% euler_hat_0 = [roll_accmag_0 + pi; -pitch_accmag_0; -yaw_accmag_0]; % TEMP!!!!!!


% q_m: Estimated-from-measurements quaternion (i.e., q from y_g (measured gyro output))
q_m = zeros(4,N_samp);
% q_m(:,1) = euler2quat(euler_0);                                                           % one of the first alternatives
% q_m(:,1) = RotMat2quat(euler2RotMat(euler_0));                                            % used till May 11, 2021; CANC?
% q_m(:,1) = RotMat2quat(euler2RotMat(euler_hat_0));                                          % used from May 11, 2021 till ...; semi-new!!!
% q_m(:,1) = RotMat2quat(rotZ(yaw_accmag_0) * rotX(pitch_accmag_0) * rotX(roll_accmag_0));    % NEW!!!!!!!!!!!!!!!!!!!!!!!!! old, canc?
q_m(:,1) = RotMat2quat(rotZ(yaw_accmag_0) * rotY(pitch_accmag_0) * rotX(roll_accmag_0));    % NEW!!!!!!!!!!!!!!!!!!!!!!!!!

for i=2:N_samp
%     % gyro bias not removed
%     q_m(:,i) = quatFirstIntegration(q_m(:,i-1), y_g(:,i), y_g(:,i-1), 'Suh');  % Possible values: 'Suh' or 'Trawny' or 'Yuan'
    
    % GYRO BIAS REMOVED!!!!!!!!
    q_m(:,i) = quatFirstIntegration(q_m(:,i-1), y_g(:,i) - b_g(:,i), y_g(:,i-1) - b_g(:,i), 'Suh');  % Possible values: 'Suh' or 'Trawny' or 'Yuan'
end

% Estimated-from-Kalman-filter quaternion initial value
q_hat__prev_prev = q_m(:,1);




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                       Some initial values (Q, x_0, P_0)
%%_________________________________________________________________________

%% Continuous-time Process noise covariance matrix, i.e., Q
% [Eq. post-10 - pre-11 Suh]
Q = [0.25 * R_g,    zeros(3),   zeros(3);
     zeros(3),      Q_b_g,      zeros(3);
     zeros(3)       zeros(3),   Q_b_a];


%% Initial state vector, i.e., x_0
q_hat_e__prev_prev = [0; 0; 0];
b_hat_g__prev_prev = [0; 0; 0];
b_hat_a__prev_prev = [0; 0; 0];
x_hat__prev_prev = [q_hat_e__prev_prev; b_hat_g__prev_prev; b_hat_a__prev_prev];


%% Initial value of the state covariance matrix, i.e., P_0
% [Eq. 5 - "How To NOT Make the Extended Kalman Filter Fail" - René Schneider]

% q_e_low = expected q_e lower bound
q_e_low = -0.02 * eye(3);                % -0.2 [rad] = -11.5 [deg]          %                                       <------------------- !
% q_e_upp = expected q_e upper bound
q_e_upp = 0.02 * eye(3);                 % 0.2 [rad] = 11.5 [deg]            %                                       <------------------- !
% q_e_0_err = q_hat_e_0 - q_e_0
q_e_0_err = (q_e_upp - q_e_low) / 2;

% b_g_low = expected b_g lower bound
b_g_low = -0.001 * eye(3);               % -0.1 [rad/s] = -5.7 [deg/s]      %                                       <------------------- !
% b_g_upp = expected b_g upper bound
b_g_upp = 0.001 * eye(3);                % 0.1 [rad/s] = 5.7 [deg/s]        %                                       <------------------- !
b_g_0_err = (b_g_upp - b_g_low) / 2;

% b_a_low = expected b_a lower bound
b_a_low = -0.02 * eye(3);                % [m / s^2]                        %                                       <------------------- !
% b_a_upp = expected b_a upper bound
b_a_upp = 0.02 * eye(3);                 % [m / s^2]                        %                                       <------------------- !
b_a_0_err = (b_a_upp - b_a_low) / 2;

% P__prev_prev = [q_e_0_err' * q_e_0_err,     zeros(3),                       zeros(3);
%                 zeros(3),                   b_g_0_err' * b_g_0_err,         zeros(3);
%                 zeros(3),                   zeros(3),                       b_a_0_err' * b_a_0_err];
% 
% 
% P__prev_prev = [0.2 * eye(3),       zeros(3),           zeros(3);
%                 zeros(3),           0.1 * eye(3),       zeros(3);
%                 zeros(3),           zeros(3),           0.1 * eye(3)];


P__prev_prev = [0.2 * eye(3),       zeros(3),           zeros(3);
                zeros(3),           0.0005 * eye(3),    zeros(3);
                zeros(3),           zeros(3),           0.05 * eye(3)];
            
% P__prev_prev = 0.05 * eye(9); % CANC!!!

%% SAVINGS
% save('DATA_interpolated___test_volo_20200213__volo_3__LOG00041_parsed_seg5.mat', 'T_theoretical', 'roll_true', 'pitch_true', 'yaw_true', 'y_g', 'y_a', 'y_m');
% save('DATA_interpolated___test_volo_20200213__volo_3__LOG00041_parsed_seg5__GLOBAL_WORKSPACE.mat');
% save('DATA_interpolated___test_volo_20200213__volo_3__LOG00041_parsed_seg5__GLOBAL_WORKSPACE_MAIN.mat');