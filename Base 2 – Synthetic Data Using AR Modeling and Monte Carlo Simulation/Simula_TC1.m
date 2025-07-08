%% ------------------------------------------------------------------------- %%
% Revised Model for Concentric Heat Exchanger                                 %
%                                                                             %
%                         GEOMETRY OF ALGETEC                                 %
%                                                                             %
% Estênio Fuzaro de Almeida                                                   %
% -------------------------------------------------------------------------- %%
%% INITIALIZATION
clearvars; close all; clc;

nfonte = 24;
marcador = 15;
num_simulations = 100; % For Monte Carlo

mySeed = 30081984;

rng_stream = RandStream('mt19937ar', 'Seed', mySeed);
RandStream.setGlobalStream(rng_stream);
%---------------------------------------------------------------------------
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual');
%--------------------------------------------------------------------------


%% Heat Exchanger Parameters
% Geometry and initial conditions of TC1
L = 1.0; % Length (m)
D_c = 28.15e-3; % Shell diameter (m)
D_t = 15.20e-3; % Tube diameter (m)
D_h = sqrt(D_c^2 - D_t^2); % Shell hydraulic diameter (m)
Area = pi * D_t * L; % Heat transfer area (m^2)
V_c = pi * (D_h / 2)^2 * L; % Fluid volume in the shell (m^3)
V_t = pi * (D_t / 2)^2 * L; % Fluid volume in the tube (m^3)

% Inlet conditions
T_t_in = 90; % Inlet temperature in the tube (°C)
T_c_in = 20; % Inlet temperature in the shell (°C)

% Flow rate and fluid properties
m_dot_t = 6.99e-2; % Mass flow rate in the tube (kg/s)
m_dot_c = 6.00e-2; % Mass flow rate in the shell (kg/s)
cp_t = 4182; % Specific heat in the tube (J/(kg·K))
cp_c = 4206; % Specific heat in the shell (J/(kg·K))
rho_t = 965.3; % Fluid density in the tube (kg/m^3)
rho_c = 998; % Fluid density in the shell (kg/m^3)

% Thermophysical properties for U calculation
k_t = 0.675; % Thermal conductivity in the tube (W/m·K)
k_c = 0.598; % Thermal conductivity in the shell (W/m·K)
mu_t = 3.15e-4; % Dynamic viscosity in the tube (Pa·s)
mu_c = 10.02e-4; % Dynamic viscosity in the shell (Pa·s)

% Calculation of Reynolds numbers
Re_t = (4 * m_dot_t) / (pi * D_t * mu_t); % Reynolds number in the tube
Re_c = (4 * m_dot_c) / (pi * D_h * mu_c); % Reynolds number in the shell

% Calculation of Prandtl numbers
Pr_t = (mu_t * cp_t) / k_t; % Prandtl number in the tube
Pr_c = (mu_c * cp_c) / k_c; % Prandtl number in the shell

% Calculation of Nusselt number
Nu_t = 0.023 * Re_t^(4/5) * Pr_t^0.3; % Nusselt number in the tube
Nu_c = 0.023 * Re_c^(4/5) * Pr_c^0.4; % Nusselt number in the shell

% Heat transfer coefficients
h_t = Nu_t * (k_t / D_t); % Convection coefficient in the tube (W/m²·K)
h_c = Nu_c * (k_c / D_h); % Convection coefficient in the shell (W/m²·K)

% Overall heat transfer coefficient
U = ((1 / h_t) + (1 / h_c))^-1; % W/(m²·K)
Umed_limpo = U;

fprintf('Calculated overall heat transfer coefficient: %.2f W/(m²·K)\n', U);

% Problem components
a1 = -Area / (V_t * rho_t * cp_t);    % Tube
a2 =  Area / (V_c * rho_c * cp_c);    % Shell
b1 = m_dot_t / (V_t * rho_t);         % Tube
b2 = m_dot_c / (V_c * rho_c);         % Shell

%% CLEAN HEAT EXCHANGER
% Integration configuration
passo = 0.01;              % Time step
tempo_total = 50;          % Simulation time (s)
intt = 0:passo:tempo_total; % Time vector
y0 = [T_t_in, T_c_in];     % Initial conditions

% Monte Carlo Simulation
for i = 1:num_simulations
    U_random_limpo(i) = unifrnd(Umed_limpo * 0.9, Umed_limpo * 1.1);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
    [t, S] = ode45(@(t,T) funcao_contracorrente_Ufixo(t,T,a1,b1,a2,b2,T_t_in,T_c_in,U_random_limpo(i)), intt, y0, opts);
    % Extracting results
    temperatura_tubo_limpo(:,i) = S(:, 1); % Tube temperature over time
    temperatura_casco_limpo(:,i) = S(:, 2); % Shell temperature over time
end

%%
% FOULING HEAT EXCHANGER
% Integration configuration
Umed_sujo = Umed_limpo * 0.75;

% Monte Carlo Simulation
temperatura_tubo_sujo = zeros(size(temperatura_tubo_limpo));
temperatura_casco_sujo = zeros(size(temperatura_casco_limpo));
U_random_sujo = zeros(1, num_simulations);
for i = 1:num_simulations
    U_random_sujo(i) = unifrnd(Umed_sujo * 0.8, Umed_sujo * 1.2);
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [t, S] = ode45(@(t, T) funcao_contracorrente_Ufixo(t, T, a1, b1, a2, b2, T_t_in, T_c_in, U_random_sujo(i)), intt, y0, opts);
    % Extracting results
    temperatura_tubo_sujo(:, i) = S(:, 1); % Tube temperature over time
    temperatura_casco_sujo(:, i) = S(:, 2); % Shell temperature over time
end

%%
% FOULING HEAT EXCHANGER 2
Umed_sujo2 = Umed_limpo * 0.4;

% Monte Carlo Simulation
temperatura_tubo_sujo2 = zeros(size(temperatura_tubo_limpo));
temperatura_casco_sujo2 = zeros(size(temperatura_casco_limpo));
U_random_sujo2 = zeros(1, num_simulations);
for i = 1:num_simulations
    U_random_sujo2(i) = unifrnd(Umed_sujo2 * 0.6, Umed_sujo2 * 1.4);
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [t, S] = ode45(@(t, T) funcao_contracorrente_Ufixo(t, T, a1, b1, a2, b2, T_t_in, T_c_in, U_random_sujo2(i)), intt, y0, opts);
    temperatura_tubo_sujo2(:, i) = S(:, 1);
    temperatura_casco_sujo2(:, i) = S(:, 2);
end

%% Saving results
save dados_AR_TC1.mat t temperatura_casco_limpo temperatura_tubo_limpo temperatura_casco_sujo temperatura_tubo_sujo temperatura_casco_sujo2 temperatura_tubo_sujo2
