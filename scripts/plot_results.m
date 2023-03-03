%% ************************************************************************
% ***********                    PLOT RESULTS                   ***********
% ***********               Author: Livio Bisogni               ***********
% *************************************************************************
%
close all
global ext_acc_detection_time_SUH ext_acc_detection_time_SAB
%
%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                INSTRUCTIONS
%__________________________________________________________________________
%
% It generates and (potentially) saves all the figures in the .png extension.
%
saving_flag = 1;                        % if '1', images (i.e., plots) are saved (as .png files);
                                        % otherwise, they ain't saved
title_flag = 1;                         % if '1', titles are added to images;
                                        % otherwise, they ain't added
%
plot_detected_ext_acc = 0;              % if '1', plot external acceleration detection instant dots in the acceleration plot
%
images_path = '/Users/a_user_name/Documents/MATLAB/source_code/images';     % It is the folder where images will be saved.
                                                                            % *** !!! change this path if needed !!! ***
main_path = '/Users/a_user_name/Documents/MATLAB/source_code';              % It is the folder containing the source code.
                                                                            % *** !!! change this path if needed !!! ***
% % images_path = '/Users/v/Documents/MATLAB/attitude_estimation/images';
% % main_path = '/Users/v/Documents/MATLAB/attitude_estimation';

%%_________________________________________________________________________




%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                            Graphics Parameters
%%_________________________________________________________________________
%
k_0 = 1;                                % First sample to be plotted; tipically, either 1 or 2
%
kkk = 0.8;                              % font-scaling multiplier factor
font_size_title = 22 * kkk;             % 20
font_size_XYlabel = 22 * kkk;           % 18
font_size_numbers_of_axes = 11 * kkk;   % -
font_size_legend = 14 * kkk;            % 14
font_size_RMStext = 14 * kkk;           % 14
%%_________________________________________________________________________

cd(images_path)





%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               ACCELERATION
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [10 800 1000 564]);

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_a(1,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), a(1,k_0:N_samples), 'r--');
    legend({'$y_{a,x}$ (Measured)', '$a_x$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif (strcmpi(sensor_data, 'real') && plot_detected_ext_acc == 1)
    for i=1:N_samples
        if (ext_acc_detection_time_SUH(i) ~= 0)
                plot(time(i), -2, 'k.');
        end
    end
    for i=1:N_samples
        if (ext_acc_detection_time_SAB(i) ~= 0)
                plot(time(i), -4, 'g.');
        end
    end
    legend({'$y_{a,x}$ (Measured)', 'Detection of External Acceleration (Suh)', 'Detection of External Acceleration (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$a_x$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_a(2,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), a(2,k_0:N_samples), 'r--');
    legend({'$y_{a,y}$ (Measured)', '$a_y$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif (strcmpi(sensor_data, 'real') && plot_detected_ext_acc == 1)
    for i=1:N_samples
        if (ext_acc_detection_time_SUH(i) ~= 0)
                plot(time(i), -2, 'k.');
        end
    end
    for i=1:N_samples
        if (ext_acc_detection_time_SAB(i) ~= 0)
                plot(time(i), -4, 'g.');
        end
    end
    legend({'$y_{a,y}$ (Measured)', 'Detection of External Acceleration (Suh)', 'Detection of External Acceleration (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$a_y$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_a(3,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), a(3,k_0:N_samples), 'r--');
    legend({'$y_{a,z}$ (Measured)', '$a_z$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif (strcmpi(sensor_data, 'real') && plot_detected_ext_acc == 1)
    for i=1:N_samples
        if (ext_acc_detection_time_SUH(i) ~= 0)
                plot(time(i), -2, 'k.');
        end
    end
    for i=1:N_samples
        if (ext_acc_detection_time_SAB(i) ~= 0)
                plot(time(i), -4, 'g.');
        end
    end
    legend({'$y_{a,z}$ (Measured)', 'Detection of External Acceleration (Suh)', 'Detection of External Acceleration (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$a_z$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)


% hAx(3) = subplot(3,1,3);
% ax = gca;
% ax.XAxis.FontSize = font_size_numbers_of_axes;
% ax.YAxis.FontSize = font_size_numbers_of_axes;
% grid on
% hold on
% plot(time(k_0:N_samples), vecnorm(y_a(:,k_0:N_samples)), 'c');
% if strcmpi(sensor_data, 'synthetic')
%     plot(time(k_0:N_samples), a(3,k_0:N_samples), 'r--');
%     legend({'$y_{a,z}$ (Measured)', '$a_z$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
% elseif (strcmpi(sensor_data, 'real') && plot_detected_ext_acc == 1)
%     for i=1:N_samples
%         if (ext_acc_detection_time_SUH(i) ~= 0)
%                 plot(time(i), -2, 'k.');
%         end
%     end
%     for i=1:N_samples
%         if (ext_acc_detection_time_SAB(i) ~= 0)
%                 plot(time(i), -4, 'g.');
%         end
%     end
%     legend({'$y_{a,z}$ (Measured)', 'Detection of External Acceleration (Suh)', 'Detection of External Acceleration (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
% end
% xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
% ylabel({'$a_z$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)


if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Acceleration', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Acceleration.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                           ACCELERATION ERRORS
%%_________________________________________________________________________
if strcmpi(sensor_data, 'synthetic')
    hAx = gobjects(1*3,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [10 800 1000 564]);
    
    hAx(1) = subplot(3,1,1);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_a(1,k_0:N_samples) - a(1,k_0:N_samples), 'c');
    % plot(time(k_0:N_samples), a_hat__next_next(1,k_0:N_samples) - a(1,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    % legend({'$\overline{a}_x - a_x$ (Measured $-$ True)', '$\widehat{a}_x - a_x$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    legend({'$y_{a,x} - a_x$ (Measured $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_x$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(2) = subplot(3,1,2);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_a(2,k_0:N_samples) - a(2,k_0:N_samples), 'c');
    % plot(time(k_0:N_samples), a_hat__next_next(2,k_0:N_samples) - a(2,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    % legend({'$\overline{a}_y - a_y$ (Measured $-$ True)', '$\widehat{a}_y - a_y$ (Kalman $-$ True)', '$Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    legend({'$y_{a,y} - a_y$ (Measured $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_y$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(3) = subplot(3,1,3);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_a(3,k_0:N_samples) - a(3,k_0:N_samples), 'c');
    % plot(time(k_0:N_samples), a_hat__next_next(3,k_0:N_samples) - a(3,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    % legend({'$\overline{a}_z - a_z$ (Measured $-$ True)', '$\widehat{a}_z - a_z$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    legend({'$y_{a,z} - a_z$ (Measured $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_z$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('Acceleration Error', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'Acceleration Error.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                       ACCELEROMETER BIAS Estimation
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [10 800 1000 564]);

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), b_hat_a__next_next(1,k_0:N_samples), 'b');
plot(time(k_0:N_samples), b_a(1,k_0:N_samples), 'r--');
legend({'$\widehat{b}_{a,x}$ (Kalman)', '$b_{a,x}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$b_{a,x}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), b_hat_a__next_next(2,k_0:N_samples), 'b');
plot(time(k_0:N_samples), b_a(2,k_0:N_samples), 'r--');
legend({'$\widehat{b}_{a,y}$ (Kalman)', '$b_{a,y}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$b_{a,y}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), b_hat_a__next_next(3,k_0:N_samples), 'b');
plot(time(k_0:N_samples), b_a(3,k_0:N_samples), 'r--');
legend({'$\widehat{b}_{a,z}$ (Kalman)', '$b_{a,z}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$b_{a,z}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Accelerometer Bias Estimation', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Accelerometer Bias Estimation.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                           EXTERNAL ACCELERATION
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [10 800 1000 564]);

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), a_b_unfiltered(1,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), a_b(1,k_0:N_samples), 'r--');
    legend({'$a_{b,x}$ (Pre-Filtering)', '$a_{b,x}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
%     plot(time(k_0:N_samples), a_b(1,k_0:N_samples), 'r--');
%     legend({'$a_{b,x}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$a_{b,x}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), a_b_unfiltered(2,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), a_b(2,k_0:N_samples), 'r--');
    legend({'$a_{b,y}$ (Pre-Filtering)', '$a_{b,y}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
%     plot(time(k_0:N_samples), a_b(2,k_0:N_samples), 'r--');
%     legend({'$a_{b,y}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$a_{b,y}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), a_b_unfiltered(3,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), a_b(3,k_0:N_samples), 'r--');
    legend({'$a_{b,z}$ (Pre-Filtering)', '$a_{b,z}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
%     plot(time(k_0:N_samples), a_b(3,k_0:N_samples), 'r--');
%     legend({'$a_{b,z}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$a_{b,z}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('External Acceleration', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'External Acceleration.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                       EXTERNAL ACCELERATION Estimation
%%_________________________________________________________________________
if strcmpi(sensor_data, 'synthetic')
    hAx = gobjects(1*3,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [10 800 1000 564]);
    
    hAx(1) = subplot(3,1,1);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), a_b(1,k_0:N_samples), 'r--');
    plot(time(k_0:N_samples), y_a(1,k_0:N_samples) - a_hat__next_next(1,k_0:N_samples) - b_hat_a__next_next(1,k_0:N_samples), 'b');           %  - b_hat_a__next_next(1,k_0:N_samples)
    if (plot_detected_ext_acc == 1)
        for i=1:N_samples
            if (ext_acc_detection_time_SUH(i) ~= 0)
                plot(time(i), -2, 'k.');
            end
        end
        for i=1:N_samples
            if (ext_acc_detection_time_SAB(i) ~= 0)
                plot(time(i), -4, 'g.');
            end
        end
        legend({'$a_{b,x}$ (True)', '$\widehat{a}_{b,x} = y_{a,x} - [1, 0, 0]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,x}$ (Kalman)', 'Detection of External Acceleration (Suh)', 'Detection of External Acceleration (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
    else
        legend({'$a_{b,x}$ (True)', '$\widehat{a}_{b,x} = y_{a,x} - [1, 0, 0]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,x}$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
    end
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{b,x}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(2) = subplot(3,1,2);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), a_b(2,k_0:N_samples), 'r--');
    plot(time(k_0:N_samples), y_a(2,k_0:N_samples) - a_hat__next_next(2,k_0:N_samples) - b_hat_a__next_next(2,k_0:N_samples), 'b');           %  - b_hat_a__next_next(2,k_0:N_samples)
    if (plot_detected_ext_acc == 1)
        for i=1:N_samples
            if (ext_acc_detection_time_SUH(i) ~= 0)
                plot(time(i), -2, 'k.');
            end
        end
        for i=1:N_samples
            if (ext_acc_detection_time_SAB(i) ~= 0)
                plot(time(i), -4, 'g.');
            end
        end
        legend({'$a_{b,y}$ (True)', '$\widehat{a}_{b,y} = y_{a,y} - [0, 1, 0]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,y}$ (Kalman)', 'Detection of External Acceleration (Suh)', 'Detection of External Acceleration (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
    else
        legend({'$a_{b,y}$ (True)', '$\widehat{a}_{b,y} = y_{a,y} - [0, 1, 0]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,y}$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
    end
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{b,y}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(3) = subplot(3,1,3);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), a_b(3,k_0:N_samples), 'r--');
    plot(time(k_0:N_samples), y_a(3,k_0:N_samples) - a_hat__next_next(3,k_0:N_samples) - b_hat_a__next_next(3,k_0:N_samples), 'b');           %  - b_hat_a__next_next(3,k_0:N_samples)
    if (plot_detected_ext_acc == 1)
        for i=1:N_samples
            if (ext_acc_detection_time_SUH(i) ~= 0)
                plot(time(i), -2, 'k.');
            end
        end
        for i=1:N_samples
            if (ext_acc_detection_time_SAB(i) ~= 0)
                plot(time(i), -4, 'g.');
            end
        end
        legend({'$a_{b,z}$ (True)', '$\widehat{a}_{b,z} = y_{a,z} - [0, 0, 1]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,z}$ (Kalman)', 'Detection of External Acceleration (Suh)', 'Detection of External Acceleration (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
    else
        legend({'$a_{b,z}$ (True)', '$\widehat{a}_{b,z} = y_{a,z} - [0, 0, 1]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,z}$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
    end
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{b,z}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('External Acceleration Estimation', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'External Acceleration Estimation.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                           EXTERNAL ACCELERATION ERROR
%%_________________________________________________________________________
if strcmpi(sensor_data, 'synthetic')
    hAx = gobjects(1*3,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [10 800 1000 564]);
    
    hAx(1) = subplot(3,1,1);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), a_b(1,k_0:N_samples), 'r--');
    % plot(time(k_0:N_samples), y_a(1,k_0:N_samples) - a_hat__next_next(1,k_0:N_samples), 'b');           %  - b_hat_a__next_next(1,k_0:N_samples)
    plot(time(k_0:N_samples), y_a(1,k_0:N_samples) - a_hat__next_next(1,k_0:N_samples) - b_hat_a__next_next(1,k_0:N_samples) - a_b(1,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$(y_{a,x} - [1, 0, 0]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,x}) - a_{b,x}$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{b,x}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(2) = subplot(3,1,2);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), a_b(2,k_0:N_samples), 'r--');
    % plot(time(k_0:N_samples), y_a(2,k_0:N_samples) - a_hat__next_next(2,k_0:N_samples), 'b');           %  - b_hat_a__next_next(2,k_0:N_samples)
    % legend({'$a_{b,y}$ (True)', '$\widehat{a}_{b,y} = y_{a,y} - [0, 1, 0]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,y}$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
    plot(time(k_0:N_samples), y_a(2,k_0:N_samples) - a_hat__next_next(2,k_0:N_samples) - b_hat_a__next_next(2,k_0:N_samples) - a_b(2,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$(y_{a,y} - [0, 1, 0]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,y}) - a_{b,y}$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{b,y}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(3) = subplot(3,1,3);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), a_b(3,k_0:N_samples), 'r--');
    % plot(time(k_0:N_samples), y_a(3,k_0:N_samples) - a_hat__next_next(3,k_0:N_samples), 'b');           %  - b_hat_a__next_next(3,k_0:N_samples)
    % legend({'$a_{b,z}$ (True)', '$\widehat{a}_{b,z} = y_{a,z} - [0, 0, 1]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,z}$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
    plot(time(k_0:N_samples), y_a(3,k_0:N_samples) - a_hat__next_next(3,k_0:N_samples) - b_hat_a__next_next(3,k_0:N_samples) - a_b(3,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$(y_{a,z} - [0, 0, 1]C(\widehat{q})\widetilde{g} - \widehat{b}_{a,z}) - a_{b,z}$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{b,z}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('External Acceleration Error', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'External Acceleration Error.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                       INTERNAL ACCELERATION Estimation
%%_________________________________________________________________________
if strcmpi(sensor_data, 'synthetic')
    hAx = gobjects(1*3,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [10 800 1000 564]);
    
    hAx(1) = subplot(3,1,1);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), y_a(1,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), a_hat__next_next(1,k_0:N_samples) + b_hat_a__next_next(1,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), a_i(1,k_0:N_samples), 'r--');
    legend({'$\widehat{a}_{i,x}$ (Kalman)', '$a_{i,x}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{i,x}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(2) = subplot(3,1,2);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), y_a(2,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), a_hat__next_next(2,k_0:N_samples) + b_hat_a__next_next(2,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), a_i(2,k_0:N_samples), 'r--');
    legend({'$\widehat{a}_{i,y}$ (Kalman)', '$a_{i,y}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{i,y}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(3) = subplot(3,1,3);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), y_a(3,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), a_hat__next_next(3,k_0:N_samples) + b_hat_a__next_next(3,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), a_i(3,k_0:N_samples), 'r--');
    legend({'$\widehat{a}_{i,z}$ (Kalman)', '$a_{i,z}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{i,z}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('Internal Acceleration Estimation', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'Internal Acceleration Estimation.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                        INTERNAL ACCELERATION ERROR
%%_________________________________________________________________________
if strcmpi(sensor_data, 'synthetic')
    hAx = gobjects(1*3,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [10 800 1000 564]);
    
    hAx(1) = subplot(3,1,1);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), y_a(1,k_0:N_samples) - a(1,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), a_hat__next_next(1,k_0:N_samples) + b_hat_a__next_next(1,k_0:N_samples) - a_i(1,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$\widehat{a}_{i,x} - a_{i,x}$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{i,x}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(2) = subplot(3,1,2);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), y_a(2,k_0:N_samples) - a(2,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), a_hat__next_next(2,k_0:N_samples) + b_hat_a__next_next(2,k_0:N_samples) - a_i(2,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$\widehat{a}_{i,y} - a_{i,y}$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{i,y}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(3) = subplot(3,1,3);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    % plot(time(k_0:N_samples), y_a(3,k_0:N_samples) - a(3,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), a_hat__next_next(3,k_0:N_samples) + b_hat_a__next_next(3,k_0:N_samples) - a_i(3,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$\widehat{a}_{i,z} - a_{i,z}$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$a_{i,z}$ [$m/s^2$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('Internal Acceleration Error', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'Internal Acceleration Error.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                             ANGULAR VELOCITY
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [10 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_g(1,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), w(1,k_0:N_samples), 'r--');
    legend({'$y_{g,x}$ (Measured)', '$\omega_x$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{g,x}$ (Measured)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\omega_x$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_g(2,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), w(2,k_0:N_samples), 'r--');
    legend({'$y_{g,y}$ (Measured)', '$\omega_y$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{g,y}$ (Measured)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\omega_y$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_g(3,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), w(3,k_0:N_samples), 'r--');
    legend({'$y_{g,z}$ (Measured)', '$\omega_z$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{g,z}$ (Measured)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\omega_z$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Angular Velocity', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Angular Velocity.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                         ANGULAR VELOCITY Estimation
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [10 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_g(1,k_0:N_samples), 'c');
plot(time(k_0:N_samples), w_hat__next_next(1,k_0:N_samples), 'b');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), w(1,k_0:N_samples), 'r--');
    legend({'$y_{g,x}$ (Measured)', '$\widehat{\omega}_x$ (Kalman)', '$\omega_x$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{g,x}$ (Measured)', '$\widehat{\omega}_x$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\omega_x$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_g(2,k_0:N_samples), 'c');
plot(time(k_0:N_samples), w_hat__next_next(2,k_0:N_samples), 'b');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), w(2,k_0:N_samples), 'r--');
    legend({'$y_{g,y}$ (Measured)', '$\widehat{\omega}_y$ (Kalman)', '$\omega_y$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{g,y}$ (Measured)', '$\widehat{\omega}_y$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\omega_y$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_g(3,k_0:N_samples), 'c');
plot(time(k_0:N_samples), w_hat__next_next(3,k_0:N_samples), 'b');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), w(3,k_0:N_samples), 'r--');
    legend({'$y_{g,z}$ (Measured)', '$\widehat{\omega}_z$ (Kalman)', '$\omega_z$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{g,z}$ (Measured)', '$\widehat{\omega}_z$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\omega_z$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Angular Velocity Estimation', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Angular Velocity Estimation.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                          ANGULAR VELOCITY ERROR
%%_________________________________________________________________________
if strcmpi(sensor_data, 'synthetic')
    hAx = gobjects(1*3,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [10 10 1000 564])
    
    hAx(1) = subplot(3,1,1);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_g(1,k_0:N_samples) - w(1,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), w_hat__next_next(1,k_0:N_samples) - w(1,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$y_{g,x} - \omega_x$ (Measured $-$ True)', '$\widehat{\omega}_x - \omega_x$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$\omega_x$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(2) = subplot(3,1,2);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_g(2,k_0:N_samples) - w(2,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), w_hat__next_next(2,k_0:N_samples) - w(2,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$y_{g,y} - \omega_y$ (Measured $-$ True)', '$\widehat{\omega}_y - \omega_y$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$\omega_y$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(3) = subplot(3,1,3);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_g(3,k_0:N_samples) - w(3,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), w_hat__next_next(3,k_0:N_samples) - w(3,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$y_{g,z} - \omega_z$ (Measured $-$ True)', '$\widehat{\omega}_z - \omega_z$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$\omega_z$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('Angular Velocity Error', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'Angular Velocity Error.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                          GYROSCOPE BIAS Estimation
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [10 800 1000 564]);

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), b_hat_g__next_next(1,k_0:N_samples), 'b');
plot(time(k_0:N_samples), b_g(1,k_0:N_samples), 'r--');
legend({'$\widehat{b}_{g,x}$ (Kalman)', '$b_{g,x}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$b_{g,x}$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), b_hat_g__next_next(2,k_0:N_samples), 'b');
plot(time(k_0:N_samples), b_g(2,k_0:N_samples), 'r--');
legend({'$\widehat{b}_{g,y}$ (Kalman)', '$b_{g,y}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$b_{g,y}$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), b_hat_g__next_next(3,k_0:N_samples), 'b');
plot(time(k_0:N_samples), b_g(3,k_0:N_samples), 'r--');
legend({'$\widehat{b}_{g,z}$ (Kalman)', '$b_{g,z}$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$b_{g,z}$ [$rad/s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Gyroscope Bias Estimation', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Gyroscope Bias Estimation.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               EULER ANGLES
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [800 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), rad2deg(roll_meas(k_0:N_samples)), 'c');
plot(time(k_0:N_samples), rad2deg(roll_hat(k_0:N_samples)), 'b');
plot(time(k_0:N_samples), rad2deg(roll_true(k_0:N_samples)), 'r--');
% plot(time(k_0:N_samples), rad2deg(roll_accmag(k_0:N_samples)), 'k-.');            % NEW (FROM ACCELEROMETER)
% legend({'$\overline{\phi}$ (Measured)', '$\widehat{\phi}$ (Kalman)', '$\phi$ (True)', '$\phi$ (from Accelerometer)'},'Interpreter','latex', 'FontSize', font_size_legend)
legend({'$\overline{\phi}$ (Measured)', '$\widehat{\phi}$ (Kalman)', '$\phi$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Roll [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), rad2deg(pitch_meas(k_0:N_samples)), 'c');
plot(time(k_0:N_samples), rad2deg(pitch_hat(k_0:N_samples)), 'b');
plot(time(k_0:N_samples), rad2deg(pitch_true(k_0:N_samples)), 'r--');
% plot(time(k_0:N_samples), rad2deg(pitch_accmag(k_0:N_samples)), 'k-.');            % NEW (FROM ACCELEROMETER)
% legend({'$\overline{\theta}$ (Measured)', '$\widehat{\theta}$ (Kalman)', '$\theta$ (True)', '$\theta$ (from Accelerometer)'},'Interpreter','latex', 'FontSize', font_size_legend)
legend({'$\overline{\theta}$ (Measured)', '$\widehat{\theta}$ (Kalman)', '$\theta$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Pitch [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), rad2deg(yaw_meas(k_0:N_samples)), 'c');
plot(time(k_0:N_samples), rad2deg(yaw_hat(k_0:N_samples)), 'b');
plot(time(k_0:N_samples), rad2deg(yaw_true(k_0:N_samples)), 'r--');
% plot(time(k_0:N_samples), rad2deg(yaw_accmag(k_0:N_samples)), 'k-.');            % NEW (FROM MAGNETOMETER)
% legend({'$\overline{\psi}$ (Measured)', '$\widehat{\psi}$ (Kalman)', '$\psi$ (True)', '$\psi$ (from Magnetometer)'},'Interpreter','latex', 'FontSize', font_size_legend)
legend({'$\overline{\psi}$ (Measured)', '$\widehat{\psi}$ (Kalman)', '$\psi$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Yaw [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Euler Angles Estimation', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Euler Angles Estimation.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               EULER ANGLES
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [800 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
% plot(time(k_0:N_samples), rad2deg(roll_meas(k_0:N_samples)), 'c');
plot(time(k_0:N_samples), rad2deg(roll_hat(k_0:N_samples)), 'b');
plot(time(k_0:N_samples), rad2deg(roll_true(k_0:N_samples)), 'r--');
% plot(time(k_0:N_samples), rad2deg(roll_accmag(k_0:N_samples)), 'k-.');            % NEW (FROM ACCELEROMETER)
% legend({'$\overline{\phi}$ (Measured)', '$\widehat{\phi}$ (Kalman)', '$\phi$ (True)', '$\phi$ (from Accelerometer)'},'Interpreter','latex', 'FontSize', font_size_legend)
legend({'$\widehat{\phi}$ (Kalman)', '$\phi$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Roll [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
% plot(time(k_0:N_samples), rad2deg(pitch_meas(k_0:N_samples)), 'c');
plot(time(k_0:N_samples), rad2deg(pitch_hat(k_0:N_samples)), 'b');
plot(time(k_0:N_samples), rad2deg(pitch_true(k_0:N_samples)), 'r--');
% plot(time(k_0:N_samples), rad2deg(pitch_accmag(k_0:N_samples)), 'k-.');            % NEW (FROM ACCELEROMETER)
% legend({'$\overline{\theta}$ (Measured)', '$\widehat{\theta}$ (Kalman)', '$\theta$ (True)', '$\theta$ (from Accelerometer)'},'Interpreter','latex', 'FontSize', font_size_legend)
legend({'$\widehat{\theta}$ (Kalman)', '$\theta$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Pitch [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
% plot(time(k_0:N_samples), rad2deg(yaw_meas(k_0:N_samples)), 'c');
plot(time(k_0:N_samples), rad2deg(yaw_hat(k_0:N_samples)), 'b');
plot(time(k_0:N_samples), rad2deg(yaw_true(k_0:N_samples)), 'r--');
% plot(time(k_0:N_samples), rad2deg(yaw_accmag(k_0:N_samples)), 'k-.');            % NEW (FROM MAGNETOMETER)
% legend({'$\overline{\psi}$ (Measured)', '$\widehat{\psi}$ (Kalman)', '$\psi$ (True)', '$\psi$ (from Magnetometer)'},'Interpreter','latex', 'FontSize', font_size_legend)
legend({'$\widehat{\psi}$ (Kalman)', '$\psi$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Yaw [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Euler Angles Estimation just KF', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Euler Angles Estimation just KF.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                         EULER ANGLES ERROR, GLOBAL
%                       NB:   wrapping is necessary
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [800 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), rad2deg(wrapToPi(roll_meas(k_0:N_samples) - roll_true(k_0:N_samples))), 'c');     % NEW (with wrapping)
plot(time(k_0:N_samples), rad2deg(wrapToPi(roll_hat(k_0:N_samples) - roll_true(k_0:N_samples))), 'b');      % NEW (with wrapping)
% plot(time(k_0:N_samples), rad2deg(roll_meas(k_0:N_samples)) - rad2deg(roll_true(k_0:N_samples)), 'c');    % without wrapping
% plot(time(k_0:N_samples), rad2deg(roll_hat(k_0:N_samples)) - rad2deg(roll_true(k_0:N_samples)), 'b');     % without wrapping
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\overline{\phi} - \phi$ (Measured $-$ True)', '$\widehat{\phi} - \phi$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Roll [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylim=get(gca,'ylim');
txt_roll_meas = ['RMSE$[\phi] =$ ', num2str(rms_roll_meas), '$^\circ$ (Measurements)'];
text(3, 8 / 10 * ylim(2), txt_roll_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
txt_roll_hat = ['RMSE$[\phi] =$ ', num2str(rms_roll_hat), '$^\circ$ (Kalman)'];
text(3, 8 / 10 * ylim(1), txt_roll_hat, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), rad2deg(wrapToPi2(pitch_meas(k_0:N_samples) - pitch_true(k_0:N_samples))), 'c');  % NEW (with wrapping)
plot(time(k_0:N_samples), rad2deg(wrapToPi2(pitch_hat(k_0:N_samples) - pitch_true(k_0:N_samples))), 'b');   % NEW (with wrapping)
% plot(time(k_0:N_samples), rad2deg(pitch_meas(k_0:N_samples)) - rad2deg(pitch_true(k_0:N_samples)), 'c');  % without wrapping
% plot(time(k_0:N_samples), rad2deg(pitch_hat(k_0:N_samples)) - rad2deg(pitch_true(k_0:N_samples)), 'b');   % without wrapping
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\overline{\theta} - \theta$ (Measured $-$ True)', '$\widehat{\theta} - \theta$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Pitch [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylim=get(gca,'ylim');
txt_pitch_meas = ['RMSE$[\theta] =$ ', num2str(rms_pitch_meas), '$^\circ$ (Measurements)'];
text(3, 8 / 10 * ylim(2), txt_pitch_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
txt_pitch_hat = ['RMSE$[\theta] =$ ', num2str(rms_pitch_hat), '$^\circ$ (Kalman)'];
text(3, 8 / 10 * ylim(1), txt_pitch_hat, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), rad2deg(wrapToPi(yaw_meas(k_0:N_samples) - yaw_true(k_0:N_samples))), 'c');       % NEW (with wrapping)
plot(time(k_0:N_samples), rad2deg(wrapToPi(yaw_hat(k_0:N_samples) - yaw_true(k_0:N_samples))), 'b');        % NEW (with wrapping)
% plot(time(k_0:N_samples), rad2deg(yaw_meas(k_0:N_samples)) - rad2deg(yaw_true(k_0:N_samples)), 'c');      % without wrapping
% plot(time(k_0:N_samples), rad2deg(yaw_hat(k_0:N_samples)) - rad2deg(yaw_true(k_0:N_samples)), 'b');       % without wrapping
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\overline{\psi} - \psi$ (Measured $-$ True)', '$\widehat{\psi} - \psi$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Yaw [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylim=get(gca,'ylim');
txt_yaw_meas = ['RMSE$[\psi] =$ ', num2str(rms_yaw_meas), '$^\circ$ (Measurements)'];
text(3, 8 / 10 * ylim(2), txt_yaw_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
txt_yaw_hat = ['RMSE$[\psi] =$ ', num2str(rms_yaw_hat), '$^\circ$ (Kalman)'];
text(3, 8 / 10 * ylim(1), txt_yaw_hat, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Euler Angles Error Global', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Euler Angles Error Global.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                       EULER ANGLES ERROR, JUST KALMAN
%                       NB:   wrapping is necessary
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [800 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
% plot(time(k_0:N_samples), rad2deg(wrapToPi(roll_meas(k_0:N_samples) - roll_true(k_0:N_samples))), 'c');     % NEW (with wrapping)
plot(time(k_0:N_samples), rad2deg(wrapToPi(roll_hat(k_0:N_samples) - roll_true(k_0:N_samples))), 'b');      % NEW (with wrapping)
% plot(time(k_0:N_samples), rad2deg(roll_meas(k_0:N_samples)) - rad2deg(roll_true(k_0:N_samples)), 'c');    % without wrapping
% plot(time(k_0:N_samples), rad2deg(roll_hat(k_0:N_samples)) - rad2deg(roll_true(k_0:N_samples)), 'b');     % without wrapping
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{\phi} - \phi$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Roll [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylim=get(gca,'ylim');
% txt_roll_meas = ['RMSE$[\phi] =$ ', num2str(rms_roll_meas), '$^\circ$ (Measurements)'];
% text(3, 8 / 10 * ylim(2), txt_roll_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
txt_roll_hat = ['RMSE$[\phi] =$ ', num2str(rms_roll_hat), '$^\circ$ (Kalman)'];
text(3, 8 / 10 * ylim(1), txt_roll_hat, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
% plot(time(k_0:N_samples), rad2deg(wrapToPi2(pitch_meas(k_0:N_samples) - pitch_true(k_0:N_samples))), 'c');  % NEW (with wrapping)
plot(time(k_0:N_samples), rad2deg(wrapToPi2(pitch_hat(k_0:N_samples) - pitch_true(k_0:N_samples))), 'b');   % NEW (with wrapping)
% plot(time(k_0:N_samples), rad2deg(pitch_meas(k_0:N_samples)) - rad2deg(pitch_true(k_0:N_samples)), 'c');  % without wrapping
% plot(time(k_0:N_samples), rad2deg(pitch_hat(k_0:N_samples)) - rad2deg(pitch_true(k_0:N_samples)), 'b');   % without wrapping
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{\theta} - \theta$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Pitch [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylim=get(gca,'ylim');
% txt_pitch_meas = ['RMSE$[\theta] =$ ', num2str(rms_pitch_meas), '$^\circ$ (Measurements)'];
% text(3, 8 / 10 * ylim(2), txt_pitch_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
txt_pitch_hat = ['RMSE$[\theta] =$ ', num2str(rms_pitch_hat), '$^\circ$ (Kalman)'];
text(3, 8 / 10 * ylim(1), txt_pitch_hat, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
% plot(time(k_0:N_samples), rad2deg(wrapToPi(yaw_meas(k_0:N_samples) - yaw_true(k_0:N_samples))), 'c');       % NEW (with wrapping)
plot(time(k_0:N_samples), rad2deg(wrapToPi(yaw_hat(k_0:N_samples) - yaw_true(k_0:N_samples))), 'b');        % NEW (with wrapping)
% plot(time(k_0:N_samples), rad2deg(yaw_meas(k_0:N_samples)) - rad2deg(yaw_true(k_0:N_samples)), 'c');      % without wrapping
% plot(time(k_0:N_samples), rad2deg(yaw_hat(k_0:N_samples)) - rad2deg(yaw_true(k_0:N_samples)), 'b');       % without wrapping
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{\psi} - \psi$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'Yaw [$^\circ$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylim=get(gca,'ylim');
% txt_yaw_meas = ['RMSE$[\psi] =$ ', num2str(rms_yaw_meas), '$^\circ$ (Measurements)'];
% text(3, 8 / 10 * ylim(2), txt_yaw_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
txt_yaw_hat = ['RMSE$[\psi] =$ ', num2str(rms_yaw_hat), '$^\circ$ (Kalman)'];
text(3, 8 / 10 * ylim(1), txt_yaw_hat, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Euler Angles Error KF', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Euler Angles Error KF.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                              MAGNETIC FIELD
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [900 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_m(1,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), m(1,k_0:N_samples), 'r--');
    legend({'$y_{m,x}$ (Measured)', '$m_x$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{m,x}$ (Measured)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$B_x [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_m(2,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), m(2,k_0:N_samples), 'r--');
    legend({'$y_{m,y}$ (Measured)', '$m_y$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{m,y}$ (Measured)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$B_y [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_m(3,k_0:N_samples), 'c');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), m(3,k_0:N_samples), 'r--');
    legend({'$y_{m,z}$ (Measured)', '$m_z$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
    legend({'$y_{m,z}$ (Measured)'},'Interpreter','latex', 'FontSize', font_size_legend)
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$B_z  [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Magnetic Field', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Magnetic Field.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                          MAGNETIC FIELD Estimation
%%_________________________________________________________________________
hAx = gobjects(1*3,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [900 10 1000 564])

hAx(1) = subplot(3,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_m(1,k_0:N_samples), 'c');
plot(time(k_0:N_samples), m_hat__next_next(1,k_0:N_samples), 'b');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), m(1,k_0:N_samples), 'r--');
    legend({'$y_{m,x}$ (Measured)', '$\widehat{m}_x$ (Kalman)', '$m_x$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
   legend({'$y_{m,x}$ (Measured)', '$\widehat{m}_x$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend) 
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$B_x [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(3,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_m(2,k_0:N_samples), 'c');
plot(time(k_0:N_samples), m_hat__next_next(2,k_0:N_samples), 'b');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), m(2,k_0:N_samples), 'r--');
    legend({'$y_{m,y}$ (Measured)', '$\widehat{m}_y$ (Kalman)', '$m_y$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
   legend({'$y_{m,y}$ (Measured)', '$\widehat{m}_y$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend) 
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$B_y [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(3,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), y_m(3,k_0:N_samples), 'c');
plot(time(k_0:N_samples), m_hat__next_next(3,k_0:N_samples), 'b');
if strcmpi(sensor_data, 'synthetic')
    plot(time(k_0:N_samples), m(3,k_0:N_samples), 'r--');
    legend({'$y_{m,z}$ (Measured)', '$\widehat{m}_z$ (Kalman)', '$m_z$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
elseif strcmpi(sensor_data, 'real')
   legend({'$y_{m,z}$ (Measured)', '$\widehat{m}_z$ (Kalman)'},'Interpreter','latex', 'FontSize', font_size_legend) 
end
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$B_z  [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Magnetic Field Estimation', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Magnetic Field Estimation.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                         MAGNETIC FIELD ERROR
%%_________________________________________________________________________
if strcmpi(sensor_data, 'synthetic')
    hAx = gobjects(1*3,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [900 10 1000 564])
    
    hAx(1) = subplot(3,1,1);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_m(1,k_0:N_samples)-m(1,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), m_hat__next_next(1,k_0:N_samples)-m(1,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$y_{m,x} - m_x$ (Measured $-$ True)', '$\widehat{m}_x - m_x$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$B_x [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(2) = subplot(3,1,2);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_m(2,k_0:N_samples)-m(2,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), m_hat__next_next(2,k_0:N_samples)-m(2,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$y_{m,y} - m_y$ (Measured $-$ True)', '$\widehat{m}_y - m_y$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$B_y [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    hAx(3) = subplot(3,1,3);
    ax = gca;
    ax.XAxis.FontSize = font_size_numbers_of_axes;
    ax.YAxis.FontSize = font_size_numbers_of_axes;
    grid on
    hold on
    plot(time(k_0:N_samples), y_m(3,k_0:N_samples)-m(3,k_0:N_samples), 'c');
    plot(time(k_0:N_samples), m_hat__next_next(3,k_0:N_samples)-m(3,k_0:N_samples), 'b');
    plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
    legend({'$y_{m,z} - m_z$ (Measured $-$ True)', '$\widehat{m}_z - m_z$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$B_z  [\mu T]$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('Magnetic Field Error', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'Magnetic Field Error.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                               QUATERNION
%%_________________________________________________________________________
hAx = gobjects(1*4,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [900 10 1000 564])

hAx(1) = subplot(4,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(1,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(1,k_0:N_samples), 'b');
plot(time(k_0:N_samples), q(1,k_0:N_samples), 'r--');
legend({'$\overline{q}_0$ (Measured)', '$\widehat{q}_0$ (Estimated)', '$q_0$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_0$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(4,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(2,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(2,k_0:N_samples), 'b');
plot(time(k_0:N_samples), q(2,k_0:N_samples), 'r--');
legend({'$\overline{q}_1$ (Measured)', '$\widehat{q}_1$ (Estimated)', '$q_1$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_1$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(4,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(3,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(3,k_0:N_samples), 'b');
plot(time(k_0:N_samples), q(3,k_0:N_samples), 'r--');
legend({'$\overline{q}_2$ (Measured)', '$\widehat{q}_2$ (Estimated)', '$q_2$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_2$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(4) = subplot(4,1,4);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(4,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(4,k_0:N_samples), 'b');
plot(time(k_0:N_samples), q(4,k_0:N_samples), 'r--');
legend({'$\overline{q}_3$ (Measured)', '$\widehat{q}_3$ (Estimated)', '$q_3$ (True)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_3$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% %%% Plot quaternion norm:
% hAx(4) = subplot(4,1,4);
% grid on
% hold on
% plot(time(k_0:N_samples), vecnorm(q_m(:,k_0:N_samples)), 'c');
% plot(time(k_0:N_samples), vecnorm(q_hat__next_next(:,k_0:N_samples)), 'b');
% plot(time(k_0:N_samples), vecnorm(q(:,k_0:N_samples)), 'r--');
% legend({'$|\overline{q}|$ (Measured)', '$|\widehat{q}|$ (Estimated)', '$|q|$ (True)'},'Interpreter','latex')
% xlabel('$t$ [$s$]')
% ylabel({'$q []$'},'Interpreter','latex')

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Quaternion Estimation', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Quaternion Estimation.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                            QUATERNION ERROR
%%_________________________________________________________________________
hAx = gobjects(1*4,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [900 10 1000 564])

hAx(1) = subplot(4,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(1,k_0:N_samples) - q(1,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(1,k_0:N_samples) - q(1,k_0:N_samples), 'b');
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\overline{q}_0 - q_0$ (Measured $-$ True)', '$\widehat{q}_0 - q_0$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_0$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(2) = subplot(4,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(2,k_0:N_samples) - q(2,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(2,k_0:N_samples) - q(2,k_0:N_samples), 'b');
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\overline{q}_1 - q_1$ (Measured $-$ True)', '$\widehat{q}_1 - q_1$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_1$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(3) = subplot(4,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(3,k_0:N_samples) - q(3,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(3,k_0:N_samples) - q(3,k_0:N_samples), 'b');
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\overline{q}_2 - q_2$ (Measured $-$ True)', '$\widehat{q}_2 - q_2$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_2$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

hAx(4) = subplot(4,1,4);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), q_m(4,k_0:N_samples) - q(4,k_0:N_samples), 'c');
plot(time(k_0:N_samples), q_hat__next_next(4,k_0:N_samples) - q(4,k_0:N_samples), 'b');
plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\overline{q}_3 - q_3$ (Measured $-$ True)', '$\widehat{q}_3 - q_3$ (Kalman $-$ True)', 'Reference'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$q_3$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% %%% Plot quaternion norm:
% hAx(4) = subplot(4,1,4);
% grid on
% hold on
% plot(time(k_0:N_samples), vecnorm(q_m(:,k_0:N_samples)) - vecnorm(q(:,k_0:N_samples)), 'c');
% plot(time(k_0:N_samples), vecnorm(q_hat__next_next(:,k_0:N_samples)) - vecnorm(q(:,k_0:N_samples)), 'b');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
% legend({'$|\overline{q}| - |q|$ (Measured $-$ True)', '$|\widehat{q}| - |q|$ (Estimated $-$ True)'},'Interpreter','latex')
% xlabel('$t$ [$s$]')
% ylabel({'$q []$'},'Interpreter','latex')

if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Quaternion Error', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Quaternion Error.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                  Process Noise Covariance Matrices (Q)
%%_________________________________________________________________________
hAx = gobjects(1*4,1);        % preallocate for the axes handles
figure('Renderer', 'painters', 'Position', [900 10 1000 564])

Q_hat_a_b_SAB_row1 = Q_hat_a_b_SAB(1,:,:);
Q_hat_a_b_SAB_row2 = Q_hat_a_b_SAB(2,:,:);
Q_hat_a_b_SAB_row3 = Q_hat_a_b_SAB(3,:,:);
Q_hat_a_b_SAB_row1 = squeeze(Q_hat_a_b_SAB_row1);
Q_hat_a_b_SAB_row2 = squeeze(Q_hat_a_b_SAB_row2);
Q_hat_a_b_SAB_row3 = squeeze(Q_hat_a_b_SAB_row3);

Q_hat_a_b_SUH_row1 = Q_hat_a_b_SUH(1,:,:);
Q_hat_a_b_SUH_row2 = Q_hat_a_b_SUH(2,:,:);
Q_hat_a_b_SUH_row3 = Q_hat_a_b_SUH(3,:,:);
Q_hat_a_b_SUH_row1 = squeeze(Q_hat_a_b_SUH_row1);
Q_hat_a_b_SUH_row2 = squeeze(Q_hat_a_b_SUH_row2);
Q_hat_a_b_SUH_row3 = squeeze(Q_hat_a_b_SUH_row3);

% 1,1
hAx(1) = subplot(9,1,1);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row1(1,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row1(1,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,1,1}}$ (Suh)', '$\widehat{Q}_{a_{b,1,1}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,1,1}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% 1,2
hAx(1) = subplot(9,1,2);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row1(2,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row1(2,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,1,2}}$ (Suh)', '$\widehat{Q}_{a_{b,1,2}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,1,2}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% 1,3
hAx(1) = subplot(9,1,3);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row1(3,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row1(3,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,1,3}}$ (Suh)', '$\widehat{Q}_{a_{b,1,3}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,1,3}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)


% 2,1
hAx(1) = subplot(9,1,4);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row2(1,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row2(1,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,2,1}}$ (Suh)', '$\widehat{Q}_{a_{b,2,1}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,2,1}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% 2,2
hAx(1) = subplot(9,1,5);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row2(2,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row2(2,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,2,2}}$ (Suh)', '$\widehat{Q}_{a_{b,2,2}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,2,2}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% 2,3
hAx(1) = subplot(9,1,6);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row2(3,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row2(3,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,2,3}}$ (Suh)', '$\widehat{Q}_{a_{b,2,3}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,2,3}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)


% 3,1
hAx(1) = subplot(9,1,7);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row3(1,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row3(1,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,3,1}}$ (Suh)', '$\widehat{Q}_{a_{b,3,1}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,3,1}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% 3,2
hAx(1) = subplot(9,1,8);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row3(2,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row3(2,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,3,2}}$ (Suh)', '$\widehat{Q}_{a_{b,3,2}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,3,2}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)

% 3,3
hAx(1) = subplot(9,1,9);
ax = gca;
ax.XAxis.FontSize = font_size_numbers_of_axes;
ax.YAxis.FontSize = font_size_numbers_of_axes;
grid on
hold on
plot(time(k_0:N_samples), Q_hat_a_b_SUH_row3(3,k_0:N_samples), 'b');
plot(time(k_0:N_samples), Q_hat_a_b_SAB_row3(3,k_0:N_samples), 'r');
% plot(time(k_0:N_samples), zeros(1,N_samples - (k_0 - 1)), 'r--');
legend({'$\widehat{Q}_{a_{b,3,3}}$ (Suh)', '$\widehat{Q}_{a_{b,3,3}}$ (Sabatini)'},'Interpreter','latex', 'FontSize', font_size_legend)
xlabel('$t$ [$s$]', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
ylabel({'$\widehat{Q}_{a_{b,3,3}}$'},'Interpreter','latex', 'FontSize', font_size_XYlabel)


if (title_flag == 1)
    % TITLE
    ax = axes;
    t1 = title('Process Noise Covariance Matrix (Q)', 'FontSize', font_size_title);
    ax.Visible = 'off'; % set(ax,'Visible','off');
    t1.Visible = 'on';  % set(t1,'Visible','on');
end
if (saving_flag == 1)
    exportgraphics(gcf, 'Q Matrices.png', 'BackgroundColor', 'white')      % crop + white background
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                                dt differences
%%_________________________________________________________________________
if strcmpi(sensor_data, 'real')
    hAx = gobjects(1*4,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [900 10 1000 564])
    
    plot(k_0:N_samples, dt_diff_real - dt);
    grid on
    hold on
    counter_dt_outliers = sum(abs(dt_diff_real - dt) > 0.0052);
    legend({'$(t_{real}(k) - t_{real}(k-1)) - T_s, k=1,...,N_{measurements}-1$'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('Sample index ($k$)', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$dt$ [$s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylim=get(gca,'ylim');
    txt_roll_meas = ['$\mathbf{card}(\{k \in 1,...,N_{measurements}-1 : |t_{real}(k) - t_{real}(k-1) - T_s| > 0.0052$', ' s\})', ' $ = $ ', num2str(counter_dt_outliers)];
    text(1200, 9 / 10 * ylim(2), txt_roll_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('dt differences', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'dt.png', 'BackgroundColor', 'white')      % crop + white background
    end
end



%% ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%                            dt cumulative error
%%_________________________________________________________________________
if strcmpi(sensor_data, 'real')
    hAx = gobjects(1*4,1);        % preallocate for the axes handles
    figure('Renderer', 'painters', 'Position', [900 10 1000 564])
    
    plot(k_0:N_samples, dt_diff_real - dt);
    grid on
    hold on
    plot(k_0:N_samples, dt_diff);
    counter_dt_outliers = sum(abs(dt_diff_real - dt) > 0.0052);
    legend({'$(t_{real}(k) - t_{real}(k-1)) - T_s, k=1,...,N_{measurements}-1$', '$t_{theoretical}(k) - t_{real}(k), k=0,...,N_{measurements}-1$'},'Interpreter','latex', 'FontSize', font_size_legend)
    xlabel('Sample index ($k$)', 'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylabel({'$dt$ [$s$]'},'Interpreter','latex', 'FontSize', font_size_XYlabel)
    ylim=get(gca,'ylim');
    txt_roll_meas = ['$\mathbf{card}(\{k \in 1,...,N_{measurements}-1 : |t_{real}(k) - t_{real}(k-1) - T_s| > 0.0052$', ' s\})', ' $ = $ ', num2str(counter_dt_outliers)];
    text(1200, 9 / 10 * ylim(2), txt_roll_meas, 'Interpreter', 'Latex', 'FontSize', font_size_RMStext)
    
    if (title_flag == 1)
        % TITLE
        ax = axes;
        t1 = title('dt cumulative error', 'FontSize', font_size_title);
        ax.Visible = 'off'; % set(ax,'Visible','off');
        t1.Visible = 'on';  % set(t1,'Visible','on');
    end
    if (saving_flag == 1)
        exportgraphics(gcf, 'dt_cumulative_error.png', 'BackgroundColor', 'white')      % crop + white background
    end
end




%%
cd(main_path)

% save('workspace_synth_200Hz_Suh_multconst_20210606.mat');