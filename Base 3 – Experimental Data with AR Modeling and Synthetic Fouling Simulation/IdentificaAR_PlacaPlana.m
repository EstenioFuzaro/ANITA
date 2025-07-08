%--------------------------------------------------------------------------
%   Identify AR Model for a Heat Exchanger - Shell Tube (Source)
%
% Author: Estênio
% Date: March 3rd 2025
%--------------------------------------------------------------------------

clc
clear
close all
nfonte = 15;

%--------------------------------------------------------------------------- 
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual');

%--------------------------------------------------------------------------
% Define stochastic parameters
rng_stream = RandStream('mt19937ar','Seed',16021971);
RandStream.setGlobalStream(rng_stream); 

%--------------------------------------------------------------------------
% Load heat exchanger data
load PlacaPlana_Data_MC.mat

% Adding noise to simulated data
fator = 0.002; % Noise intensity factor
t = Time_Vector;

n_exp = 5;
n_sim = 45;
n_fou = 50;

% Separate datasets
temperatura_fria_limpo_exp = All_Curves(:, 1:n_exp)';  % Experimental clean curves
temperatura_fria_limpo_sim = All_Curves(:, n_exp+1:n_exp+n_sim)'; % Simulated clean curves
temperatura_fria_sujo = All_Curves(:, n_exp+n_sim+1:end)';    % Fouled curves

% Compute standard deviation for correct noise scaling
std_clean = std(temperatura_fria_limpo_sim, 0, 2); % Compute std per row
std_fouled = std(temperatura_fria_sujo, 0, 2); 

% Generate noise with the same shape as the data
noise_clean = fator * std_clean .* randn(size(temperatura_fria_limpo_sim));
noise_fouled = fator * std_fouled .* randn(size(temperatura_fria_sujo));

% Add noise to the simulated data
temperatura_fria_limpo_sim = temperatura_fria_limpo_sim + noise_clean;
temperatura_fria_sujo = temperatura_fria_sujo + noise_fouled;

% Combine clean experimental and simulated data
temperatura_fria_limpo = [temperatura_fria_limpo_exp; temperatura_fria_limpo_sim];

%--------------------------------------------------------------------------
% AR Model Identification using AIC

dt = t(2) - t(1); % Sampling rate

% Select a random clean data sample for training
np = randperm(n_exp+n_sim);
np_train = np(1);  
y = temperatura_fria_limpo(np_train, :);  % Training data

% Normalize the training data (zero mean and unit variance)
y = zscore(y);

% Format for MATLAB System Identification Toolbox
data = iddata(y', [], dt); % Convert to IDDATA format

%--------------------------------------------------------------------------
% Choosing the model order A(q)y(k) = e(k) using Akaike Information Criterion

na = 2:2:20;  % Orders of the AR model to evaluate
Akaike = zeros(size(na)); % Store AIC values

for i = 1:length(na)
    model = ar(data, na(i));
    Akaike(i) = aic(model); % Compute Akaike Information Criterion
end

%--------------------------------------------------------------------------
% Plot Akaike Criterion Results
figure()
stem(na,Akaike)
xlabel('Order')
ylabel('AIC')

figure()
stem(na, Akaike, 'filled', 'k');
hold on
plot(na, Akaike, '-r', 'LineWidth', 1.5); % Smooth connection
xlabel('AR Model Order', 'Interpreter', 'latex', 'FontSize', nfonte);
ylabel('AIC Value', 'Interpreter', 'latex', 'FontSize', nfonte);
title('Akaike Criterion for AR Model Order Selection', 'Interpreter', 'latex', 'FontSize', nfonte);
grid on;
set(gca, 'FontSize', nfonte, 'TickLabelInterpreter', 'latex');

% Save Akaike data for future analysis
%save('akaike_data_he1.mat', 'na', 'Akaike');

%--------------------------------------------------------------------------
% Reference Model
ordem = 6;
model_ref = ar(data, ordem);
coeficiente_limpo = polydata(model_ref); % Extracting AR model coefficients

%--------------------------------------------------------------------------
% Validation Stage
% Select a **different** sample for validation
np_val = np(2);

x_v = temperatura_fria_limpo(np_val, :); % Validation data

% Normalize validation data
x_v = zscore(x_v);

% Format validation data
data_v = iddata(x_v', [], dt);

%--------------------------------------------------------------------------
% Model Validation - Compare AR Model with Experimental Data
figure()
compare(model_ref, data_v, 1)

% Extract and replot comparison results
lines = findobj(gcf, 'Type', 'line');
xData = get(lines, 'XData');
yData = get(lines, 'YData');

figure();
hold on;
for i = length(xData):-1:2
    if i == 3
        plot(xData{i}, yData{i}, 'b', 'LineWidth', 1.5); % Experimental data in blue
    elseif i == 2
        plot(xData{i}, yData{i}, 'r', 'LineWidth', 1.5); % Model AR data in red
    else
        plot(xData{i}, yData{i}); % Other data
    end
end
hold off;
box on;
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', nfonte);
ylabel('Shell Temperature [$^\circ$C]', 'Interpreter', 'latex', 'FontSize', nfonte);
xlim([0 30]);
legend('Experimental', 'AR Model', 'Location', 'Best');

%--------------------------------------------------------------------------
% Residual Analysis
figure()
resid(data_v, model_ref)

% Compute the residuals
residuals = resid(data_v, model_ref);

% Compute the auto-correlation function (ACF) of residuals
[acf, lags] = autocorr(residuals.OutputData, 'NumLags', 20);

% Compute the 99% confidence bounds
n = length(residuals.OutputData);
confBounds = 2.576 / sqrt(n); % 99% confidence level

% Plot ACF with confidence bounds
figure;
plot(lags, acf, 'b', 'LineWidth', 1.5); % ACF plot in blue
hold on;
plot(lags, confBounds * ones(size(lags)), 'r--', 'LineWidth', 1.5); % Upper bound
plot(lags, -confBounds * ones(size(lags)), 'r--', 'LineWidth', 1.5); % Lower bound
hold off;
xlabel('Lags', 'Interpreter', 'latex', 'FontSize', nfonte);
ylabel('Autocorrelation', 'Interpreter', 'latex', 'FontSize', nfonte);
legend('ACF', '99% Confidence Interval');

%--------------------------------------------------------------------------
% One-step-ahead Prediction Error
erro = pe(data_v, model_ref);
erro_ref = erro.OutputData;

figure()
plot(t, erro_ref, 'k', 'LineWidth', 1.5);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', nfonte);
ylabel('Prediction Error', 'Interpreter', 'latex', 'FontSize', nfonte);
title('One-Step-Ahead Prediction Error', 'Interpreter', 'latex', 'FontSize', nfonte);
grid on;

% Save the validated model and data
%save model_clean_data_source model_ref data_v

%--------------------------------------------------------------------------
% Calculating Prediction Errors for Unknown Conditions (clean or fouled)
%--------------------------------------------------------------------------
error = [];

% Data known to be clean
clean_indices = randperm(n_exp+n_sim-1) + 1; % Shuffle indices for clean data, starting from 2 to 50

for i = 1:n_exp+n_sim-1
   x = temperatura_fria_limpo(:, clean_indices(i)); % Use shuffled indices
   % x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x = iddata(x, [], dt);
   erro = pe(data_x, model_ref);
   erro = erro.OutputData;

   % Storing errors
   error = [error erro];
end

% Data known to be fouled - continue filling the error matrix
fouled_indices = randperm(n_exp+n_sim); % Shuffle indices

for i = 1:n_fou
   x = temperatura_fria_sujo(:, fouled_indices(i)); % Use shuffled indices
   % x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x=iddata(x,[],dt);
   erro = pe(data_x,model_ref);
   erro = erro.OutputData;

   % Storing errors
   error= [error erro];

end

%--------------------------------------------------------------------------
% Calculating Index X_1
%--------------------------------------------------------------------------
% Calculating standard deviation in reference condition
REF = var(erro_ref);

% Calculating error variance in unknown condition (clean or fouled)
UNK = var(error);

X_1 = UNK./REF;

nfonte = 18;

figure()
set(gcf, 'Position', [100, 100, 800, 400]) % Set rectangular format (width 800, height 400)
h1 = plot(X_1(1:(n_exp+n_sim-1)), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); hold on
h2 = plot(n_exp+n_sim:n_exp+n_sim+n_fou-1,X_1(n_exp+n_sim:n_exp+n_sim+n_fou-1), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % Add markers with fill
%plot(X_1(200:299), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r') % Add markers with fill
min_val = min(X_1(n_exp+n_sim-1));
yline(min_val, '--r', 'LineWidth', 2, 'DisplayName', 'Min Value (1:49)') % Red dashed line
grid on % Add grid
xlabel('Observations', 'Interpreter', 'latex', 'FontSize', nfonte) % X label with LaTeX
ylabel('Prediction Error Ratio', 'Interpreter', 'latex', 'FontSize', nfonte) % Y label with LaTeX
legend('Clean','Fouled')
set(gca, 'FontSize', nfonte, 'TickLabelInterpreter', 'latex') % Axis configuration
%save x1_data_source X_1

%--------------------------------------------------------------------------
% Calculating X_2 - ratio between variance of AR model coefficients

A1 = [];
% Data known to be clean
for i =2:(n_exp+n_sim)

   x = temperatura_fria_limpo(:,np(i));
   %x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x=iddata(x,[],dt);
   model = ar(data_x,ordem);
   coeficientes = polydata(model);

   % Storing errors
   A1= [A1; coeficientes];

end

ordem_sujo = ordem+2;
% Data known to be fouled
A2 = [];
for i =1:n_fou

   x = temperatura_fria_sujo(:,i);
   %x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x=iddata(x,[],dt);
   model = ar(data_x,ordem_sujo);
   coeficientes = polydata(model);

   % Storing errors
   A2= [A2; coeficientes];

end

% Validate and check one of the fouled models
x = temperatura_fria_sujo(:,np(1));
%x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
x = zscore(x);
data_x=iddata(x,[],dt);

figure()
compare(model,data_x,1)

figure()
resid(model,data_x)
%save model_fouled_data_source model data_x

%--------------------------------------------------------------------------
% Visualizing X_1 by X_2

% Calculation of standard deviation in reference condition
REF = var(coeficiente_limpo);

% Calculation of error variance in unknown condition (clean or fouled)
UNK1 = var(A1');
UNK2 = var(A2');

X_2 = [UNK1 UNK2]./REF;

figure()
set(gcf, 'Position', [100, 100, 800, 600]) % Set rectangular format (width 800, height 600)
% Plot of the first 99 points
plot(X_1(1:(n_exp+n_sim-1)), X_2(1:(n_exp+n_sim-1)), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b') % Light blue filled
hold on
% Plot of points 100 to 199
plot(X_1((n_exp+n_sim):(n_exp+n_sim+n_fou-1)), X_2((n_exp+n_sim):(n_exp+n_sim+n_fou-1)), 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', 'r') % Red X
% Plot of the last points (200)
%plot(X_1(200:299), X_2(200:299), 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k') % Black filled pentagon
% Grid and labels
grid on % Add grid for better visualization
xlabel('Feature 1', 'Interpreter', 'latex', 'FontSize', nfonte) % X label with LaTeX
ylabel('Feature 2', 'Interpreter', 'latex', 'FontSize', nfonte) % Y label with LaTeX
set(gca, 'FontSize', nfonte, 'TickLabelInterpreter', 'latex') % Axis configuration
% Add legend
legend({'Clean', 'Fouled'}, ...
    'Interpreter', 'latex', 'FontSize', nfonte, 'Location', 'best')
%save features_data_source X_1 X_2

% % Save in EPS and JPG formats
% saveas(gcf, 'FEATURES.eps', 'epsc') % Save as EPS
% saveas(gcf, 'FEATURES.jpg') % Save as JPG

%---------------------------------------------------------------------------
% Nova feature: Energia das curvas de temperatura
%---------------------------------------------------------------------------

curva_identificacao = np(1); % Curva usada para identificar o modelo

% Remover a curva selecionada do conjunto de curvas limpas
temperatura_fria_limpo_corrigido = temperatura_fria_limpo;
temperatura_fria_limpo_corrigido(curva_identificacao,:) = [];

% Calcular a energia para cada curva com ruído
energia_limpo = trapz(t, temperatura_fria_limpo_corrigido', 1); % Integral para curvas limpas corrigidas com ruído
energia_sujo1 = trapz(t, temperatura_fria_sujo', 1);  % Integral para curvas sujas 1 com ruído

% Criar índices para as curvas
indice_limpo = 1:(n_exp+n_sim-1);
indice_sujo1 = (n_exp+n_sim):(n_exp+n_sim+n_fou-1);

% Visualizar as energias com scatter plots
figure;

% Scatter para curvas limpas (azul)
scatter(indice_limpo, energia_limpo, (n_exp+n_sim), 'b', 'filled', 'DisplayName', 'Clean');
hold on;

% Scatter para curvas sujas 1 (vermelho)
scatter(indice_sujo1, energia_sujo1, n_fou, 'r', 'filled', 'DisplayName', 'Fouled');

% Scatter para curvas sujas 2 (vermelho com borda preta)
%scatter(indice_sujo2, energia_sujo2, 50, 'r', 'filled', ...
%    'MarkerEdgeColor', 'k', 'DisplayName', 'Sujo 2');

% Configurar o gráfico
grid on;
xlabel('Observations', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Energy (Integral)', 'Interpreter', 'latex', 'FontSize', 14);
%title('Nova Feature: Energia das Curvas de Temperatura', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

% Para salvar a energia como feature X_2, descomentar:
X_2 = [energia_limpo energia_sujo1];
%
%-------------------------------------------------------------------
% Saving features X_1 and X_2 and model
%-------------------------------------------------------------------
save('features_PlacaPlana.mat', 'X_1', 'X_2'); % Save to .mat file
save('model_PlacaPlana','model_ref','coeficiente_limpo','erro_ref','ordem','ordem_sujo')
