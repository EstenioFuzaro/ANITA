%--------------------------------------------------------------------------
%   Identify AR Model for a Heat Exchanger - Target
%
% Author: Samuel (UNESP)
% Date: 01/07/2025
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
% define stochastic parameters
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setGlobalStream(rng_stream); 

%--------------------------------------------------------------------------
% Load heat exchanger data
load dados_AR_TC2

% Non-parametric analysis

figure()

% Plot the first set of data with legend
h1 = plot(t, temperatura_casco_limpo(:,1)', 'g'); hold on
h2 = plot(t, temperatura_casco_sujo(:,1)', 'm');
legend([h1, h2], 'Clean', 'Fouled') % Include only the first two handles in the legend

% Plot the remaining data without adding them to the legend
plot(t, temperatura_casco_limpo(:,2:end)', 'g', 'HandleVisibility', 'off');
plot(t, temperatura_casco_sujo(:,2:end)', 'm', 'HandleVisibility', 'off');

% Axis labels and configuration
xlabel('Time [s]', 'Interpreter', 'latex', 'FontSize', nfonte) % X label with LaTeX
ylabel('Shell Temperature [$^\circ$C]', 'Interpreter', 'latex', 'FontSize', nfonte) % Y label with LaTeX
set(gca, 'FontSize', nfonte, 'TickLabelInterpreter', 'latex') % Axis configuration


%--------------------------------------------------------------------------
% 20 shell temperature data - clean and fouled (with 25% reduction of
% overall heat transfer coefficient U)

dt = t(2)-t(1); % sampling rate

% clean data - for reference model
np = randperm(20);  % randomly selecting one of the 20 data
y = temperatura_casco_limpo(:,np(1));
% adding % RMS white noise to the signal
fator = 0.002;
y = y + fator*sqrt(mean(y.^2))*randn(length(t),1);

y = zscore(y);

% All other data will be assumed unknown

% Format for ID Toolbox 
data = iddata(y,[],dt); % reference condition

%--------------------------------------------------------------------------
% Choosing the model order A(q)y(k) = e(k)

na = 2:2:20;    % order of polynomial A(q)

for i=1:length(na)

    model = ar(data,na(i));
    Akaike(i) = aic(model);
end

figure()
stem(na,Akaike)
xlabel('Order')
ylabel('AIC')
save akaike_data_he2 na Akaike

%--------------------------------------------------------------------------
% Reference Model
ordem = 10;
model_ref = ar(data,ordem);
coeficiente_limpo = polydata(model_ref); % Extracting AR model coefficients

%--------------------------------------------------------------------------
% Validation stage
% selecting another data to validate (different from the one used to identify
% the model)
% adding % RMS noise to the signal
x_v = temperatura_casco_limpo(:,np(2));

% adding % RMS white noise to the signal
x_v = x_v + fator*sqrt(mean(x_v.^2))*randn(length(t),1);
x_v = zscore(x_v);
data_v=iddata(x_v,[],dt);

% Generate the figure using the compare function
figure()
compare(model_ref, data_v, 1)

% Find the lines in the current figure
lines = findobj(gcf, 'Type', 'line');

% Extract the data from the lines
xData = get(lines, 'XData');
yData = get(lines, 'YData');

% Plot the extracted data as you wish, starting from the second index
figure();
hold on;
for i = length(xData):-1:2
    if i == 3
        plot(xData{i}, yData{i}, 'b',LineWidth=1.5); % Plot the second valid data set in blue
    elseif i == 2
        plot(xData{i}, yData{i},'r',LineWidth=1.5); % Plot the first valid data set in red
    else
        plot(xData{i}, yData{i}); % Plot any additional data sets with default color
    end
end
hold off;
box on;
% xlabel('Tempo [s]');
% ylabel('Temperatura do Casco [$^\circ$C]');
% xlim([0 50]);
% legend('Experimental', 'Modelo AR', 'Location', 'Best'); % Adjust the legend as needed
xlabel('Time [s]');
ylabel('Shell Temperature [$^\circ$C]');
xlim([0 50]);
legend('Experimental', 'AR Model', 'Location', 'Best'); % Adjust the legend as needed

% Residual analysis
figure()
resid(data_v, model_ref)

% Extract residuals from the model
residuals = resid(data_v, model_ref);

% Calculate the auto-correlation function (ACF) of the residuals
[acf, lags] = autocorr(residuals.OutputData, 'NumLags', 20);

% Calculate the 99% confidence bounds
n = length(residuals.OutputData);
confBounds = 2.576 / sqrt(n); % 2.576 is the z-value for 99% confidence

% Plot the ACF with confidence bounds as lines
figure;
plot(lags, acf, 'b', 'LineWidth', 1.5); % Plot ACF as a blue line
hold on;
plot(lags, confBounds * ones(size(lags)), 'r--', 'LineWidth', 1.5); % Upper confidence bound as a red dashed line
plot(lags, -confBounds * ones(size(lags)), 'r--', 'LineWidth', 1.5); % Lower confidence bound as a red dashed line
hold off;
% xlabel('Atrasos');
% ylabel('Auto-Correla\c{c}\~ao');
% legend('Auto-Correla\c{c}\~ao', 'Intervalo de Confian\c{c}a (99$\%$)');
xlabel('Lags');
ylabel('Autocorrelation');
legend('Autocorrelation', 'Confidence Interval (99$\%$)');

% One-step-ahead prediction error
erro = pe(data_v,model_ref);
erro_ref = erro.OutputData;

figure()
plot(t,erro_ref)
xlabel('Time (s)');
ylabel('Prediction Error')
save model_clean_data_target.mat

%--------------------------------------------------------------------------
% Calculating Prediction Errors for Unknown Conditions (clean or fouled)
%--------------------------------------------------------------------------
error = [];

% Data known to be clean
for i =2:20

   x = temperatura_casco_limpo(:,np(i));
   x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x=iddata(x,[],dt);
   erro = pe(data_x,model_ref);
   erro = erro.OutputData;

   % Storing errors
   error= [error  erro];

end

% Data known to be fouled - continue filling the error matrix

for i = 1:20

   x = temperatura_casco_sujo(:,i);
   x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x=iddata(x,[],dt);
   erro = pe(data_x,model_ref);
   erro = erro.OutputData;

   % Storing errors
   error= [error erro];

end


for i = 1:20

   x = temperatura_casco_sujo2(:,i);
   x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
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
load model_TC1.mat erro_ref
% Calculating standard deviation in reference condition
REF = var(erro_ref);

% Calculating error variance in unknown condition (clean or fouled)
UNK = var(error);

X_1 = UNK./REF;

nfonte = 18;

figure()
set(gcf, 'Position', [100, 100, 800, 400]) % Set rectangular format (width 800, height 400)
h1 = plot(X_1(1:19), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b'); hold on
h2 = plot(20:39,X_1(20:39), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r'); % Add markers with fill
%plot(X_1(200:299), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'r') % Add markers with fill
min_val = min(X_1(1:19));
yline(min_val, '--r', 'LineWidth', 2, 'DisplayName', 'Min Value (1:99)') % Red dashed line
grid on % Add grid
xlabel('Observations', 'Interpreter', 'latex', 'FontSize', nfonte) % X label with LaTeX
ylabel('Prediction Error Ratio', 'Interpreter', 'latex', 'FontSize', nfonte) % Y label with LaTeX
legend('Clean','Fouled')
set(gca, 'FontSize', nfonte, 'TickLabelInterpreter', 'latex') % Axis configuration

%--------------------------------------------------------------------------
% Calculating X_2 - ratio between variance of AR model coefficients

A1 = [];
% Data known to be clean
for i =2:20

   x = temperatura_casco_limpo(:,np(i));
   x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
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
for i =1:20

   x = temperatura_casco_sujo(:,i);
   x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x=iddata(x,[],dt);
   model = ar(data_x,ordem_sujo);
   coeficientes = polydata(model);

   % Storing errors
   A2= [A2; coeficientes];

end

% Validate and check one of the fouled models
x = temperatura_casco_sujo(:,3);
x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
x = zscore(x);
data_x=iddata(x,[],dt);

figure()
compare(model,data_x,1)

figure()
resid(model,data_x)
save model_fouled_data_target model data_x

ordem_sujo2 = ordem_sujo+2;
A3 = [];
for i =1:20

   x = temperatura_casco_sujo2(:,i);
   x = x + fator*sqrt(mean(x.^2))*randn(length(t),1);
   %x = zscore(x);
   data_x=iddata(x,[],dt);
   model = ar(data_x,ordem_sujo2);
   coeficientes = polydata(model);

   % Storing errors
   A3= [A3; coeficientes];

end

%--------------------------------------------------------------------------
% Visualizing X_1 by X_2

load model_TC1.mat coeficiente_limpo
% Calculation of standard deviation in reference condition
REF = var(coeficiente_limpo);

% Calculation of error variance in unknown condition (clean or fouled)
UNK1 = var(A1');
UNK2 = var(A2');
UNK3 = var(A3');

X_2 = [UNK1 UNK2 UNK3]./REF;

figure()
set(gcf, 'Position', [100, 100, 800, 600]) % Set rectangular format (width 800, height 600)

% Plot of the first 19 points
plot(X_1(1:19), X_2(1:19), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b') % Light blue filled
hold on

% Plot of points 20 to 39
plot(X_1(20:39), X_2(20:39), 'x', 'MarkerSize', 10, 'LineWidth', 1.5, 'Color', 'r') % Red X

% Plot of the last points (40)
%plot(X_1(40:59), X_2(40:59), 'p', 'MarkerSize', 12, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k') % Black filled pentagon

% Grid and labels
grid on % Add grid for better visualization
xlabel('Feature 1', 'Interpreter', 'latex', 'FontSize', nfonte) % X label with LaTeX
ylabel('Feature 2', 'Interpreter', 'latex', 'FontSize', nfonte) % Y label with LaTeX
set(gca, 'FontSize', nfonte, 'TickLabelInterpreter', 'latex') % Axis configuration

% Add legend
legend({'Clean', 'Fouled'}, ...
    'Interpreter', 'latex', 'FontSize', nfonte, 'Location', 'best')

% % Save in EPS and JPG formats
% saveas(gcf, 'FEATURES.eps', 'epsc') % Save as EPS
% saveas(gcf, 'FEATURES.jpg') % Save as JPG

%---------------------------------------------------------------------------
% Nova feature: Energia das curvas de temperatura
%---------------------------------------------------------------------------

curva_identificacao = np(1); % Curva usada para identificar o modelo

% Remover a curva selecionada do conjunto de curvas limpas
temperatura_casco_limpo_corrigido = temperatura_casco_limpo;
temperatura_casco_limpo_corrigido(:, curva_identificacao) = [];

% Adicionar ruído às curvas
fator = 0.002; % Fator de ruído definido anteriormente
temperatura_casco_limpo_com_ruido = temperatura_casco_limpo_corrigido + ...
    fator * sqrt(mean(temperatura_casco_limpo_corrigido.^2, 1)) .* randn(size(temperatura_casco_limpo_corrigido));
temperatura_casco_sujo_com_ruido = temperatura_casco_sujo + ...
    fator * sqrt(mean(temperatura_casco_sujo.^2, 1)) .* randn(size(temperatura_casco_sujo));
temperatura_casco_sujo2_com_ruido = temperatura_casco_sujo2 + ...
    fator * sqrt(mean(temperatura_casco_sujo2.^2, 1)) .* randn(size(temperatura_casco_sujo2));

% Calcular a energia para cada curva com ruído
energia_limpo = trapz(t, temperatura_casco_limpo_com_ruido, 1); % Integral para curvas limpas corrigidas com ruído
energia_sujo1 = trapz(t, temperatura_casco_sujo_com_ruido, 1);  % Integral para curvas sujas 1 com ruído
energia_sujo2 = trapz(t, temperatura_casco_sujo2_com_ruido, 1); % Integral para curvas sujas 2 com ruído

% Criar índices para as curvas
indice_limpo = 1:19;
indice_sujo1 = 20:39;
indice_sujo2 = 40:59;

% Visualizar as energias com scatter plots
figure;

% Scatter para curvas limpas (azul)
scatter(indice_limpo, energia_limpo, 50, 'b', 'filled', 'DisplayName', 'Limpo');
hold on;

% Scatter para curvas sujas 1 (vermelho)
scatter(indice_sujo1, energia_sujo1, 50, 'r', 'filled', 'DisplayName', 'Sujo 1');

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
X_2 = [energia_limpo energia_sujo1 energia_sujo2];
%
%-------------------------------------------------------------------
% Saving features X_1 and X_2 and model
%-------------------------------------------------------------------
save('features_TC2.mat', 'X_1', 'X_2'); % Save to .mat file
save('model_TC2','model_ref','coeficiente_limpo','erro_ref','ordem','ordem_sujo','ordem_sujo2')
