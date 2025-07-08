% Clear workspace
clear; close all; clc;

% Define parameters
num_realizations = 5;  % Number of experimental curves
num_simulations = 45;  % Number of synthetic normal heat exchanger curves
num_fouled = 50;       % Number of synthetic fouled heat exchanger curves
t_max = 30;            % Maximum simulation time (seconds)
num_time_points = 100; % Fixed number of time points for consistency

% Variance Scaling Factor (Lower values reduce variance)
variance_scaling_factor = 0.3; % Adjust between 0.3 and 0.7 as needed

% Define time vector (uniform for all curves)
Time_Vector = linspace(0, t_max, num_time_points)';

% Preallocate storage matrices
Exp_Curves = NaN(num_time_points, num_realizations);   % Experimental curves
Sim_Curves = NaN(num_time_points, num_simulations);    % Normal heat exchanger curves
Fouled_Curves = NaN(num_time_points, num_fouled);      % Fouled heat exchanger curves
All_Curves = NaN(num_time_points, num_realizations + num_simulations + num_fouled); % Combined dataset

% Initialize parameter storage
k_vals = zeros(num_realizations, 1);
tau_vals = zeros(num_realizations, 1);
T_initial_vals = zeros(num_realizations, 1);

figure; hold on;

% Process experimental data
for i = 1:num_realizations
    % Load data
    filename = sprintf('PlacaPlana_E%d.mat', i);
    data = load(filename);
    
    % Extract time and temperature
    t_i = data.t;
    T_i = data.Temp_s_fria;

    % Normalize time to start at zero
    t_i = t_i - t_i(1);

    % Trim data at 30 seconds
    valid_idx = t_i <= t_max;
    t_i = t_i(valid_idx);
    T_i = T_i(valid_idx);

    % Interpolate onto fixed time grid
    T_interp = interp1(t_i, T_i, Time_Vector, 'linear', 'extrap');

    % Store data
    Exp_Curves(:, i) = T_interp;
    All_Curves(:, i) = T_interp;
    T_initial_vals(i) = T_interp(1);

    % Estimate steady-state gain and time constant
    k_est = mean(T_interp(end-5:end));
    model_fun = @(p, t) p(1) * (1 - exp(-t / p(2))) + p(3);
    p0 = [k_est, max(t_i)/3, T_interp(1)];
    opts = optimset('Display', 'off');
    p_fit = lsqcurvefit(model_fun, p0, Time_Vector, T_interp, [0, 0, -Inf], [Inf, Inf, Inf], opts);

    k_vals(i) = p_fit(1);
    tau_vals(i) = p_fit(2);

    % Plot experimental curves (blue)
    plot(Time_Vector, T_interp, 'b', 'LineWidth', 1.5);
end

% Compute mean and scaled standard deviation
k_mean = mean(k_vals);
k_std_adj = variance_scaling_factor * std(k_vals);
tau_mean = mean(tau_vals);
tau_std_adj = variance_scaling_factor * std(tau_vals);
T_initial_mean = mean(T_initial_vals);
T_initial_std = std(T_initial_vals);

% Estimate Gamma distribution parameters
alpha_k = (k_mean^2) / (k_std_adj^2);
theta_k = (k_std_adj^2) / k_mean;
alpha_tau = (tau_mean^2) / (tau_std_adj^2);
theta_tau = (tau_std_adj^2) / tau_mean;
alpha_Ti = (T_initial_mean^2) / (T_initial_std^2);
theta_Ti = (T_initial_std^2) / T_initial_mean;

%% Monte Carlo Simulation (Normal Heat Exchanger)
for j = 1:num_simulations
    % Sample parameters from Gamma distribution (with controlled variance)
    k_sim = gamrnd(alpha_k, theta_k);
    tau_sim = gamrnd(alpha_tau, theta_tau);
    T_initial_sim = gamrnd(alpha_Ti, theta_Ti);

    % Generate temperature response
    T_sim = k_sim * (1 - exp(-Time_Vector / tau_sim)) + T_initial_sim;

    % Store and plot simulated curves
    Sim_Curves(:, j) = T_sim;
    All_Curves(:, num_realizations + j) = T_sim;
    plot(Time_Vector, T_sim, 'g', 'LineWidth', 0.5);
end

%% Generate Fouled Heat Exchanger Curves
% Adjusted mean values for fouled heat exchanger
k_fouled_mean = 0.7 * k_mean;
tau_fouled_mean = 1.25 * tau_mean;

% Adjust standard deviations for fouled heat exchanger
k_std_fouled_adj = variance_scaling_factor * std(k_vals);
tau_std_fouled_adj = variance_scaling_factor * std(tau_vals);

% Compute adjusted Gamma parameters
alpha_k_fouled = (k_fouled_mean^2) / (k_std_fouled_adj^2);
theta_k_fouled = (k_std_fouled_adj^2) / k_fouled_mean;
alpha_tau_fouled = (tau_fouled_mean^2) / (tau_std_fouled_adj^2);
theta_tau_fouled = (tau_std_fouled_adj^2) / tau_fouled_mean;

% Generate fouled heat exchanger curves
for j = 1:num_fouled
    k_fouled = gamrnd(alpha_k_fouled, theta_k_fouled);
    tau_fouled = gamrnd(alpha_tau_fouled, theta_tau_fouled);
    T_initial_fouled = gamrnd(alpha_Ti, theta_Ti);

    % Simulate temperature response
    T_fouled = k_fouled * (1 - exp(-Time_Vector / tau_fouled)) + T_initial_fouled;

    % Store and plot
    Fouled_Curves(:, j) = T_fouled;
    All_Curves(:, num_realizations + num_simulations + j) = T_fouled;
    plot(Time_Vector, T_fouled, 'r', 'LineWidth', 0.5);
end

%% Save Data
save('PlacaPlana_Data_MC.mat', 'Exp_Curves', 'Sim_Curves', 'Fouled_Curves', 'All_Curves', 'Time_Vector');

disp('Data saved successfully to PlacaPlana_Data_MC.mat');

%% Plot Everything Together
figure; hold on;

% Plot normal simulated heat exchanger curves (Green, background)
for j = 1:num_simulations
    plot(Time_Vector, Sim_Curves(:, j), 'g', 'LineWidth', 0.5);
end

% Plot fouled heat exchanger curves (Red, middle layer)
for j = 1:num_fouled
    plot(Time_Vector, Fouled_Curves(:, j), 'r', 'LineWidth', 0.5);
end

% Plot experimental curves (Blue, on top)
for i = 1:num_realizations
    plot(Time_Vector, Exp_Curves(:, i), 'b', 'LineWidth', 1.5);
end

xlabel('Time (s)');
ylabel('Temperature');
title('Experimental, Normal, and Fouled Heat Exchanger Curves');
grid on;
legend({'Simulated (Normal)', 'Simulated (Fouled)', 'Experimental'}, 'Location', 'best');
