%% ************************************************************************
% ***********              GENERATE SYNTHETIC DATA              ***********
% ***********               Author: Livio Bisogni               ***********
% *************************************************************************
%
%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                INSTRUCTIONS
%__________________________________________________________________________
%
% It generates the synthetic data.
%
% Please select one of the following rate_mode:
%                   * 'null':       null rate
%                   * 'const':      constant rate
%                   * 'roll':       zero(0-50 s) - phi=phi+90°     - zero (150-200 s)
%                   * 'pitch':      zero(0-50 s) - theta=theta+45° - zero (150-200 s)
%                   * 'yaw':        zero(0-50 s) - psi=psi+90°     - zero (150-200 s)
%                   * 'mult_const': miscellaneous (multiple) constant rates
%                   * 'mult_ramp':  miscellaneous (multiple), fast, ramp-like rates
%
rate_mode = 'mult_ramp';
%
% Please choose either 'ON' or 'OFF' (case insensitive, btw):
%
noise_gyro = 'on';  % if 'ON', it adds noise to gyroscope measurements
noise_acc  = 'on';  % if 'ON', it adds noise to accelerometer measurements
noise_mag  = 'on';  % if 'ON', it adds noise to magnetometer measurements
%
bias_gyro  = 'on';  % if 'ON', it adds bias to gyroscope measurements
bias_acc   = 'on';  % if 'ON', it adds bias to accelerometer measurements
%
ext_acc    = 'on';  % if 'ON', it creates 3 external accelerations:
                    %   * [10; 5; 20] from 80 to 81 s
                    %   * [0; -7; 0] from 120 to 122 s
                    %   * [-4; -3; 8] from 140 to 140.5 s
%%_________________________________________________________________________




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                              GLOBAL VARIABLES
%%_________________________________________________________________________
global dt N_samples         % time-related stuff
global w y_g                % true (theoretical) and measured (bias and noise added) gyroscope data, respectively
global a y_a                % true (theoretical) and measured (bias and noise added) accelerometer data, respectively
global m y_m                % true (theoretical) and measured (bias and noise added) magnetometer data, respectively

global R_a R_m                      % covariance measurement matrices of the accelerometer / magnetometer
global g_standard g_tilde m_tilde
global M_1 M_2 gamma
global Q                            % Continuous covariance process noise matrix
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
%                                 CONSTANTS
%%_________________________________________________________________________
% [Eq. post-3 Suh]
% g_standard: Standard acceleration of gravity
g_tilde_intensity = 9.80665;  % g_standard                                         % [m / s^2]                             <-------------------
g_tilde_intensity = 9.805185;  % g at Viareggio                                         % [m / s^2]                             <-------------------
% g_tilde: Gravitational field, expressed in the navigation frame
g_tilde = [0; 0; g_tilde_intensity];

% [Eq. post-3 Suh]
% Data source:      https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml#igrfwmm
% dip_angle: Dip angle (aka magnetic inclination) @ Viareggio (43° 52' 52" N, 10° 14' 6" E), IT, February 11th, 2021
% P.S.: (43° 52' 52" N, 10° 14' 6" E) = (43.88111° N, 10.23500° E)
dip_angle = deg2rad(60.10803); % 60° 6' 29"                     % [rad]                                 <-------------------
% m_tilde_intensity: Magnetic field intensity
m_tilde_intensity = 47.1179;                                    % [microT]                              <-------------------
% m_tilde: Magnetic field (simpler model), expressed in the navigation frame
m_tilde = [cos(dip_angle); 0; -sin(dip_angle)] * m_tilde_intensity;
% m_tilde: Magnetic field (wicked model) predicted @ Viareggio (43° 52' 52" N, 10° 14' 6" E), IT, February 11th, 2021, expressed in the navigation frame
% m_tilde = [23.4485; 1.2545; -40.8497];                          % [microT]                              <-------------------


% Adaptive Suh's algorithm parameters
% [Eq. 36 Suh]
Q_b_g = 0.000001 * eye(3);                                      % 0.000001                              <------------------- !
Q_b_a = 0.000001 * eye(3);                                      % 0.000001                              <------------------- !
M_1 = 3;                                                        % []                                    <------------------- (not used)
M_2 = 2;                                                        % []                                    <------------------- (not used)
gamma = 0.1;                                                    %                                       <------------------- (not used)


% Norm-based Sabatini's algorithm parameters (used by Suh)
% [Eq. 37 Suh]
s = 10;                                                         % 10; try also 1 and 100                <-------------------
epsilon_a = 0.25;                                               % [m / s^2]                             <-------------------




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                   TIME
%%_________________________________________________________________________

% f_sampling: Sampling frequency	(defined in: [post-Eq. 151 Trawny] - [pre-Eq. 152 Trawny])
%                                   NB: dt = t(k+1) - t(k)
f_sampling = 200;                                               % [Hz]                                  <------------------- !

% dt: Sampling time; it is CONSTANT
dt = 1 / f_sampling;                                            % [s]

% sim_duration: Simulation duration
sim_duration = 200;                                             % [s]                                   <-------------------

% N_samples: Number of samples
N_samples = sim_duration / dt;                                  % []
N_samples__theoretical = N_samples;                             % []

% time: Time stamp array [..., time(k-1), time(k), time(k+1), ...]
time = zeros(N_samples, 1);                                     % [s]
time(1) = 0;
for i = 2:N_samples
    time(i) = (i - 1) * dt;
end

% ext_acc_detection_time: = 1 only when external acceleration is detected; otherwise (in all other time steps) it remains null
ext_acc_detection_time_SAB = zeros(1,N_samples);
ext_acc_detection_time_SUH = zeros(1,N_samples);




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               GYROSCOPE
%%_________________________________________________________________________

% sigma_g: Standard deviation of the gyroscope global errors (noises, external perturbations, calibration defects etc.)
if strcmpi(noise_gyro, 'ON')
    sigma_g = 0.006;  % 0.01;    0.06;   0.001;                 % [rad / s]                             <------------------- !
elseif strcmpi(noise_gyro, 'OFF')
    sigma_g = 0;                                                % [rad / s]
else
    error('Invalid noise_gyro. Possible values: ''ON'' or ''OFF''.');
end

% n_g:  Gyroscope noise
n_g_mean = 0;                                                   % [rad / s]
n_g_variance = sigma_g^2;                                       % [rad^2 / s^2]
n_g = zeros(3, N_samples);                                      % [rad / s]
for i = 1:N_samples
    % [Eq. 133-134 Trawny]
    n_g(:,i) = n_g_mean + randn(3, 1) * sqrt(n_g_variance);
end

% R_g: Covariance measurement matrix of the gyroscope:
R_g = n_g_variance * eye(3);                                    % [rad^2 / s^2]


% n_b_g:  Gyro bias (b_g) driving noise; it corresponds to the process noise [random Walk process]:
n_b_g_mean = 0;                                                 % [rad / s]
sigma_b_g = 0;                                                  % [rad / s]                             <------------------- !
n_b_g_variance = sigma_b_g^2;                                   % [rad^2 / s^2]
n_b_g = zeros(3, N_samples);                                    % [rad / s]
for i = 1:N_samples
    % [Eq. 136-137 Trawny]
    n_b_g(:,i) = n_b_g_mean + randn(3, 1) * sqrt(n_b_g_variance);
end

% b: Gyro bias
b_g = zeros(3, N_samples);                                      % [rad / s]
if strcmpi(bias_gyro, 'ON')
    b_g(:,1) = [-0.019; 0.013; -0.006]; % [-0.019; 0.013; -0.006] [-0.02; 0.01; -0.06]; [-0.0067; -0.0049; -0.0097];     %                                       <------------------- !
elseif strcmpi(bias_gyro, 'OFF')
    % b_g is already void
else
    error('Invalid bias_gyro. Possible values: ''ON'' or ''OFF''.');
end
for i = 2:N_samples
    % [Eq. 135 Trawny]
    b_g(:,i) = b_g(:,i-1) + n_b_g(:,i);
end


%% w: True (theoretical) Angular Velocity
w = zeros(3, N_samples);                                        % [rad / s]

if strcmpi(rate_mode, 'null')
    %% Null rate
    for i = 1:N_samples
        w(:,i) =[0; 0; 0];
    end

elseif strcmpi(rate_mode, 'const')
    %% Constant rate
    w_x = 0.1;                                                  % [rad / s]
    w_y = 0.2;                                                  % [rad / s]
    w_z = 0.3;                                                  % [rad / s]
    
    for i = 1:N_samples
        w(:,i) =[w_x; w_y; w_z];
    end
    
elseif strcmpi(rate_mode, 'roll')
    %% Roll
    for i = 1:N_samples
        if (1 <= i && i <= 5/20 * N_samples)
            w(:,i) =[0; 0; 0];
        elseif (5/20 * N_samples + 1 <= i && i <= 15/20 * N_samples)
            w(:,i) =[deg2rad(0.9); 0; 0];
%           w(:,i) =[0.1; 0; 0];
        elseif (15/20 * N_samples + 1 <= i && i <= N_samples)
            w(:,i) =[0; 0; 0];
        end
    end

elseif strcmpi(rate_mode, 'pitch')
    %% Pitch
    for i = 1:N_samples
        if (1 <= i && i <= 5/20 * N_samples)
            w(:,i) =[0; 0; 0];
        elseif (5/20 * N_samples + 1 <= i && i <= 15/20 * N_samples)
            w(:,i) =[0; deg2rad(0.45); 0];
        elseif (15/20 * N_samples + 1 <= i && i <= N_samples)
            w(:,i) =[0; 0; 0];
        end
    end

elseif strcmpi(rate_mode, 'yaw')
    %% Yaw
    for i = 1:N_samples
        if (1 <= i && i <= 5/20 * N_samples)
            w(:,i) =[0; 0; 0];
        elseif (5/20 * N_samples + 1 <= i && i <= 15/20 * N_samples)
            w(:,i) =[0; 0; deg2rad(0.9)];
        elseif (15/20 * N_samples + 1 <= i && i <= N_samples)
            w(:,i) =[0; 0; 0];
        end
    end
    
elseif strcmpi(rate_mode, 'mult_const')
    %% Miscellaneous (multiple) constant rates
    for i = 1:N_samples
        if (1 <= i && i <= 5/20 * N_samples)
            w(:,i) =[0; 0; 0];
            
        elseif (5/20 * N_samples + 1 <= i && i <= 7/20 * N_samples)
            w(:,i) =[0.1; 0; 0];
        elseif (7/20 * N_samples + 1 <= i && i <= 8/20 * N_samples)
            w(:,i) =[0; 0; 0];
            
        elseif (8/20 * N_samples + 1 <= i && i <= 10/20 * N_samples)
            w(:,i) =[0; 0.1; 0];
        elseif (10/20 * N_samples + 1 <= i && i <= 11/20 * N_samples)
            w(:,i) =[0; 0; 0];
            
        elseif (11/20 * N_samples + 1 <= i && i <= 13/20 * N_samples)
            w(:,i) =[0; 0; 0.1];
        elseif (13/20 * N_samples + 1 <= i && i <= 14/20 * N_samples)
            w(:,i) =[0; 0; 0];
            
        elseif (14/20 * N_samples + 1 <= i && i <= 17/20 * N_samples)
            w(:,i) =[-0.3; 0.1; -0.2];
        
        elseif (17/20 * N_samples + 1 <= i && i <= N_samples)
            w(:,i) =[0; 0; 0];
        end
    end
    
elseif strcmpi(rate_mode, 'mult_ramp')
    %% Miscellaneous (multiple), fast, ramp-like rates
    for i = 1:N_samples
        if (1 <= i && i <= 5/20 * N_samples)
            w(:,i) =[0; 0; 0];
        
        elseif (5/20 * N_samples + 1 <= i && i <= 7/20 * N_samples)
            w(:,i) =[(i - 5/20 * N_samples + 1) / N_samples * 20; 0; 0];
        elseif (7/20 * N_samples + 1 <= i && i <= 8/20 * N_samples)
            w(:,i) =[0; 0; 0];
            
        elseif (8/20 * N_samples + 1 <= i && i <= 10/20 * N_samples)
            w(:,i) =[0; (i - 8/20 * N_samples + 1) / N_samples * 10; 0];
        elseif (10/20 * N_samples + 1 <= i && i <= 11/20 * N_samples)
            w(:,i) =[0; 0; 0];
            
        elseif (11/20 * N_samples + 1 <= i && i <= 13/20 * N_samples)
            w(:,i) =[0; 0; (i - 11/20 * N_samples + 1) / N_samples * 10];
        elseif (13/20 * N_samples + 1 <= i && i <= 14/20 * N_samples)
            w(:,i) =[0; 0; 0];
            
        elseif (14/20 * N_samples + 1 <= i && i <= 17/20 * N_samples)
            w(:,i) =[(i - 8/20 * N_samples + 1) / N_samples * (10); (i - 8/20 * N_samples + 1) / N_samples * (-30); (i - 8/20 * N_samples + 1) / N_samples * (-20)];
        
        elseif (17/20 * N_samples + 1 <= i && i <= N_samples)
            w(:,i) =[0; 0; 0];
        end
    end
else
    warning('Wrong value for rate_mode. Possible values: ''null'', ''const'', ''roll'', ''pitch'', ''yaw'', ''mult_const'' or ''mult_ramp''.')
end


% [Eq. 3b Suh] ([Eq. 132 Trawny])
% y_g: Measured Angular Velocity (i.e., gyro output)           % [rad / s]
% w:    True (theoretical) Angular Velocity
% b_g:  Gyroscope bias
% n_g:  Gyroscope noise
y_g = w + b_g + n_g;




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%           TRUE INITIAL EULER ANGLES & TRUE QUATERNION & TRUE C_n_b
%%_________________________________________________________________________

% Initial true Euler angles
roll_0  = -20;                                                  % [deg]                                 <------------------- !
pitch_0 = 15;                                                   % [deg]                                 <------------------- !
yaw_0   = 30;                                                   % [deg]                                 <------------------- !

% roll_0  = 0;                                                    % [deg]                                 <------------------- !
% pitch_0 = 0;                                                    % [deg]                                 <------------------- !
% yaw_0   = 0;                                                    % [deg]                                 <------------------- !

% roll_0 = rand_range(-180,180);                                      % [deg]
% pitch_0 = rand_range(-90,90);                                       % [deg]
% yaw_0 = rand_range(-180,180);                                       % [deg]

roll_0 = deg2rad(roll_0);
pitch_0 = deg2rad(pitch_0);
yaw_0 = deg2rad(yaw_0);

euler_0 = [roll_0; pitch_0; yaw_0];

% q: true quaternion, i.e., q from w (true omega)
q = zeros(4,N_samples);
% q(:,1) = euler2quat(euler_0);                                                             % one of the first alternatives
q(:,1) = RotMat2quat(euler2RotMat(euler_0));                                                % used from May 11, 2021 till ...; semi-new!!!
q(:,1) = RotMat2quat(rotZ(yaw_0) * rotY(pitch_0) * rotX(roll_0));                        	% NEW!!!!!!!!!!!!!!!!!!!!!!!!!

for i=2:N_samples
    q(:,i) = quatFirstIntegration(q(:,i-1), w(:,i), w(:,i-1), 'Suh');  % Possible values: 'Suh' or 'Trawny' or 'Yuan'
end

% [Eq. 1 Suh] ([Eq. 79 Trawny])
% C_n_b: True associated transformation matrix (function of the true attitude quaternion q)
C_n_b = zeros(3,3,N_samples);
for i = 1:N_samples
    C_n_b(:,:,i) = quat2RotMat(q(:,i));
end



    
%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               ACCELEROMETER
%%_________________________________________________________________________

% sigma_a: Standard deviation of the accelerometer global errors (noises, external perturbations, calibration defects etc.)
if strcmpi(noise_acc, 'ON')
    sigma_a = 0.048; %0.004;     0.039;                         % [m / s^2]                             <------------------- !
elseif strcmpi(noise_acc, 'OFF')
    sigma_a = 0;                                                % [m / s^2]
else
    error('Invalid noise_acc. Possible values: ''ON'' or ''OFF''.');
end

% n_a: Accelerometer noise
n_a_mean = 0;                                                   % [m / s^2]
n_a_variance = sigma_a^2;                                       % [m^2 / s^4]
n_a = zeros(3, N_samples);
for i = 1:N_samples
    % [Eq. 133-134 Trawny]
    n_a(:,i) = n_a_mean + randn(3, 1) * sqrt(n_a_variance);
end

% R_a: Covariance measurement matrix of the accelerometer:
R_a = n_a_variance * eye(3);                                    % [m^2 / s^4]


% n_b_a:  Accelerometer bias (b_a) driving noise; it corresponds to the process noise [random Walk process]:
n_b_a_mean = 0;                                                 % [m / s^2]
sigma_b_a = 0;                                                  % [m / s^2]                             <------------------- !
n_b_a_variance = sigma_b_a^2;                                   % [m^2 / s^4]
n_b_a = zeros(3, N_samples);                                    % [m / s^2]
for i = 1:N_samples
    % [Eq. 136-137 Trawny]
    n_b_a(:,i) = n_b_a_mean + randn(3, 1) * sqrt(n_b_a_variance);
end

% b_a: Accelerometer bias
b_a = zeros(3, N_samples);                                      % [m / s^2]
if strcmpi(bias_acc, 'ON')
    b_a(:,1) = [0.07; 0.033; -0.044]; % [0.07; 0.035; -0.04];  [0.09; -0.01; 0.7];                      <------------------- !
elseif strcmpi(bias_acc, 'OFF')
    % b_a is already void
else
    error('Invalid bias_acc. Possible values: ''ON'' or ''OFF''.');
end
for i = 2:N_samples
    % [Eq. 135 Trawny]
    b_a(:,i) = b_a(:,i-1) + n_b_a(:,i);
end
b_a_hat__prev_prev = b_a(:,1);


% a_b: (Unknown and unwanted) External Acceleration
%   The unwanted external acceleration information (a_b) in the accelerometer will be estimated.
%   This acceleration is caused by the change of velocity in magnitude and/or direction.
%   The external acceleration is the main source of errors in attitude estimation.
a_b = zeros(3,N_samples);                                       % [m / s^2]
%
% Generate 3 spikes of external accelerations
if strcmpi(ext_acc, 'ON')
    a_b = create_extAcc(a_b, 80, [10; 5; 20], 1);
    a_b = create_extAcc(a_b, 120, [0; -7; 0], 1);
    a_b = create_extAcc(a_b, 140, [-4; -3; 8], 1);
    
    a_b_unfiltered = a_b;
    span = 1 * f_sampling;                              % window size (number of samples)
    a_b(1,:) = smooth(a_b(1,:), span);
    a_b(2,:) = smooth(a_b(2,:), span);
    a_b(3,:) = smooth(a_b(3,:), span);

elseif strcmpi(ext_acc, 'OFF')
    % a_b is already void
else
    error('Invalid ext_acc. Possible values: ''ON'' or ''OFF''.');
end


% ****** ALMOST USELESS (just for plotting purposes) ******
% a_i: True (theoretical) internal acceleration                 % [m / s^2]
a_i = zeros(3,N_samples);
for i = 1:N_samples
    a_i(:,i) = C_n_b(:,:,i) * g_tilde;
end
% *********************************************************


% a: True (theoretical) acceleration                            % [m / s^2]
a = zeros(3,N_samples);
for i = 1:N_samples
    a(:,i) = C_n_b(:,:,i) * g_tilde + a_b(:,i);
end


% [Eq. 3b Suh]
% y_a  Measured acceleration (i.e., accelerometer output)      % [m / s^2]
% b_a:  Accelerometer bias
% n_a:  Accelerometer noise
y_a = a + b_a + n_a;




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               MAGNETOMETER
%%_________________________________________________________________________

% sigma_m: Standard deviation of the magnetometer global errors (noises, external perturbations, calibration defects etc.)
if strcmpi(noise_mag, 'ON')
    sigma_m = 2;          %0.2;    [0.022; 0.111; 0.028];       % [microT]                             <------------------- !
elseif strcmpi(noise_mag, 'OFF')
    sigma_m = 0;                                                % [microT]
else
    error('Invalid noise_mag. Possible values: ''ON'' or ''OFF''.');
end

% n_m: Magnetometer noise
n_m_mean = 0;                                                   % [microT]
n_m_variance = sigma_m^2;                                       % [(microT)^2]
n_m = zeros(3, N_samples);                                      % [microT]
for i = 1:N_samples
    % [Eq. 133-134 Trawny]
    n_m(:,i) = n_m_mean + randn(3, 1) * sqrt(n_m_variance);
end

% R_m: Covariance measurement matrix of the magnetometer:
R_m = n_m_variance * eye(3);                                    % [(microT)^2]


% m: True (theoretical) magnetic field                          % [microT]
m = zeros(3,N_samples);
for i = 1:N_samples
    m(:,i) = C_n_b(:,:,i) * m_tilde;
end


% [Eq. 3c Suh]
% y_m: Measured magnetic field (i.e., Magnetometer output)
% n_m:  Magnetometer noise
y_m = m + n_m;                                                  % [microT]




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                     INITIAL ESTIMATED EULER ANGLES AND
%              ESTIMATED QUATERNION BY INTEGRATING ANGULAR VELOCITY
%%_________________________________________________________________________

% Initial Euler Angles, estimated from Accelerometer & Magnetometer measurements:
% [https://habr.com/en/post/499190/]
%
k_0 = 2;                                                        % Index of the first non-null measurement; tipically, either 1 or 2
% Estimated initial Roll (from Accelerometer measurements)
% [Eq. 25 - Tilt Sensing Using a Three-Axis Accelerometer]
roll_accmag_0  = (atan2(y_a(2,k_0), y_a(3,k_0)));                                   % [rad]
%
% Estimated initial Pitch (from Accelerometer measurements)
% [Eq. 26 - Tilt Sensing Using a Three-Axis Accelerometer]
pitch_accmag_0 = -wrapToPi(atan2(y_a(1,k_0), sqrt(y_a(2,k_0).^2 + y_a(3,k_0).^2))); % [rad]
%
% Estimated initial Yaw (from Magnetometer measurements)
% [Eq. 22 - Implementing a Tilt-Compensated eCompass using Accelerometer and Magnetometer Sensors]
C_n_b__0 = euler2RotMat([roll_accmag_0; pitch_accmag_0; 0]);
C_n_b__0 = rotZ(0) * rotY(pitch_accmag_0) * rotX(roll_accmag_0);                      % !!!!!!!!!!!!!!!

C_b_n__0 = C_n_b__0';
y_m_0_NEW = C_b_n__0 * y_m(:,k_0);
yaw_accmag_0 = -atan2(y_m_0_NEW(2), y_m_0_NEW(1));                                  % [rad]

% Estimated initial Euler angles
euler_hat_0 = [roll_accmag_0; pitch_accmag_0; yaw_accmag_0];


% q_m: Estimated-from-measurements quaternion (i.e., q from y_g (measured gyro output))
q_m = zeros(4,N_samples);
% q_m(:,1) = euler2quat(euler_0);                                                           % one of the first alternatives
% q_m(:,1) = RotMat2quat(euler2RotMat(euler_0));                                            % used till May 11, 2021; CANC?
q_m(:,1) = RotMat2quat(euler2RotMat(euler_hat_0));                                          % used from May 11, 2021 till ...; semi-new!!!
q_m(:,1) = RotMat2quat(rotZ(yaw_accmag_0) * rotY(pitch_accmag_0) * rotX(roll_accmag_0));    % NEW!!!!!!!!!!!!!!!!!!!!!!!!!

for i=2:N_samples
    q_m(:,i) = quatFirstIntegration(q_m(:,i-1), y_g(:,i), y_g(:,i-1), 'Suh');  % Possible values: 'Suh' or 'Trawny' or 'Yuan'
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
q_e_low = -0.2 * eye(3);                % -0.2 [rad] = -11.5 [deg]          %                                       <------------------- !
% q_e_upp = expected q_e upper bound
q_e_upp = 0.2 * eye(3);                 % 0.2 [rad] = 11.5 [deg]            %                                       <------------------- !
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

P__prev_prev = [q_e_0_err' * q_e_0_err,     zeros(3),                       zeros(3);
                zeros(3),                   b_g_0_err' * b_g_0_err,         zeros(3);
                zeros(3),                   zeros(3),                       b_a_0_err' * b_a_0_err];