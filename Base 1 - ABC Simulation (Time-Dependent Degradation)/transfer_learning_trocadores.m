%% ------------------------------------------------------------------------- %%
% Trocador de calor concêntrico - contra-corrente                             %
%                                                                             %
%                    Classificação SVM e Transfer Learning                    %
%                                                                             %
% Vitória Batista Godoy e Estênio Fuzaro de Almeida                           %
% -------------------------------------------------------------------------- %%
%% INICIALIZAÇÃO
clearvars; close all; clc;

nfonte = 15;
marcador=5;

mySeed = 160403;

rng_stream = RandStream('mt19937ar', 'Seed', mySeed);
RandStream.setGlobalStream(rng_stream);
%---------------------------------------------------------------------------
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot,'defaultAxesZMinorGrid','on','defaultAxesZMinorGridMode','manual');
set(groot, 'defaultTextFontSize', nfonte);
set(groot, 'defaultAxesFontSize', nfonte);
set(groot, 'defaultLegendFontSize', nfonte);

% Load data
load('he1_data_1year.mat');
load('he1_data_4year.mat');
load('he2_data_1year.mat');
load('he2_data_4year.mat');

%% Passo 1: Organizar Dados do TC1
X_limpo_tc1 = [noisy_temp_tube_he1_1year', noisy_temp_shell_he1_1year'];
X_incrustado_tc1 = [noisy_temp_tube_he1_4year', noisy_temp_shell_he1_4year'];
X_tc1 = [X_limpo_tc1; X_incrustado_tc1];
y_tc1 = [ones(size(X_limpo_tc1, 1), 1); -ones(size(X_incrustado_tc1, 1), 1)];

%% Passo 2: Organizar Dados do TC2 (Somente Dados de Teste)
X_limpo_tc2 = [noisy_temp_tube_he2_1year', noisy_temp_shell_he2_1year'];
X_incrustado_tc2 = [noisy_temp_tube_he2_4year', noisy_temp_shell_he2_4year'];
X_tc2 = [X_limpo_tc2; X_incrustado_tc2];
y_tc2 = [ones(size(X_limpo_tc2, 1), 1); -ones(size(X_incrustado_tc2, 1), 1)];

%% Passo 3: Treinar SVM com Dados do TC1 Apenas
SVMModel_tc1 = fitcsvm(X_tc1, y_tc1, 'KernelFunction', 'linear', 'Standardize', true);

%% Passo 4: Plotar Espaço de Características Original com Fronteira de Decis\~ao
figure;
gscatter(X_limpo_tc1(:,1), X_limpo_tc1(:,2), ones(size(X_limpo_tc1, 1), 1), 'b', 'x', [], 'off');
hold on;
gscatter(X_incrustado_tc1(:,1), X_incrustado_tc1(:,2), -ones(size(X_incrustado_tc1, 1), 1), 'r', 'o', [], 'off');
gscatter(X_limpo_tc2(:,1), X_limpo_tc2(:,2), ones(size(X_limpo_tc2, 1), 1), 'g', '+', [], 'off');
gscatter(X_incrustado_tc2(:,1), X_incrustado_tc2(:,2), -ones(size(X_incrustado_tc2, 1), 1), 'm', 's', [], 'off');
svmBoundary(X_tc1, SVMModel_tc1);  % Função auxiliar para plotar a fronteira de decisão
hold off;
xlabel('Temperatura do Tubo [$^{\circ}$C]');
ylabel('Temperatura do Casco [$^{\circ}$C]');
legend({'Fonte - Limpo', 'Fonte - Incrustado', ...
        'Alvo - Limpo', 'Alvo - Incrustado', ...
        'Fronteira de Decis\~ao'}, 'Location', 'Best');

%% Passo 5: Aplicar Transferência de Aprendizado com JDA
X_tc1_normalizado = zscore(X_tc1);
X_tc2_normalizado = zscore(X_tc2);

lambda = 1;
dim = 2;
kernel_type = 'linear';
[Zs_tc1, Zt_tc2] = JDA(X_tc1_normalizado, X_tc2_normalizado, y_tc1, lambda, dim, kernel_type);

%% Passo 6: Treinar SVM com Dados Adaptados do TC1
SVMModel_jda = fitcsvm(Zs_tc1, y_tc1, 'KernelFunction', 'linear', 'Standardize', true);

% Passo 7: Plotar Espaço de Características Latentes com Fronteira de Decisão (JDA)
figure;
hold on;

% Plotar dados adaptados do TC1 (Limpo: azul 'x', Incrustado: vermelho 'o')
gscatter(Zs_tc1(y_tc1 == 1, 1), Zs_tc1(y_tc1 == 1, 2), ones(sum(y_tc1 == 1), 1), 'b', 'x', [], 'off');
gscatter(Zs_tc1(y_tc1 == -1, 1), Zs_tc1(y_tc1 == -1, 2), -ones(sum(y_tc1 == -1), 1), 'r', 'o', [], 'off');

% Plotar dados adaptados do TC2 (Limpo: verde '+', Incrustado: magenta 's')
gscatter(Zt_tc2(y_tc2 == 1, 1), Zt_tc2(y_tc2 == 1, 2), ones(sum(y_tc2 == 1), 1), 'g', '+', [], 'off');
gscatter(Zt_tc2(y_tc2 == -1, 1), Zt_tc2(y_tc2 == -1, 2), -ones(sum(y_tc2 == -1), 1), 'm', 's', [], 'off');

% Usar a mesma função de fronteira para consistência
svmBoundary([Zs_tc1; Zt_tc2], SVMModel_jda);

hold off;
xlabel('Atributo Adaptado 1');
ylabel('Atributo Adaptado 2');
legend({'Fonte - Limpo', 'Fonte - Incrustado', ...
        'Alvo - Limpo', 'Alvo - Incrustado', ...
        'Fronteira de Decis\~ao'}, 'Location', 'Best');

% Ajustar a visibilidade das linhas dos eixos
set(gca, 'Box', 'on', 'Color', 'none', 'XColor', 'k', 'YColor', 'k');  % k para preto

%% Passo 8: Gerar Previsões e Matrizes de Confusão
% Previsões sem Transferência de Aprendizado
y_pred = predict(SVMModel_tc1, [X_tc1; X_tc2]);
confMat = confusionmat([y_tc1; y_tc2], y_pred);

% Previsões com Transferência de Aprendizado (JDA)
y_pred_jda = [predict(SVMModel_jda, Zs_tc1); predict(SVMModel_jda, Zt_tc2)];
confMat_jda = confusionmat([y_tc1; y_tc2], y_pred_jda);

%% Passo 9: Definir Colormap para Matrizes de Confusão
numColors = 100;
blueColormap = [linspace(0.678, 0, numColors)', linspace(0.847, 0, numColors)', linspace(0.902, 0.545, numColors)'];

%% Passo 10: Plotar Matriz de Confusão - Antes da Transferência de Aprendizado
figure;
confMat_percentage = confMat ./ sum(confMat, 2) * 100;
imagesc(confMat_percentage);
colormap(blueColormap);
colorbar('TickLabelInterpreter', 'latex');
xlabel('Classe Prevista');
ylabel('Classe Verdadeira');
%title('Matriz de Confusão Antes da Transferência de Aprendizado');
axis square;
grid off;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Incrustado', 'Limpo'});
set(gca, 'YTick', 1:2, 'YTickLabel', {'Incrustado', 'Limpo'});

% Anotar células com valores em porcentagem
textStrings = num2str(confMat_percentage(:), '%.1f\\%%');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:2);
text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');

%% Passo 11: Plotar Matriz de Confusão - Após a Transferência de Aprendizado (JDA)
figure;
confMat_jda_percentage = confMat_jda ./ sum(confMat_jda, 2) * 100;
imagesc(confMat_jda_percentage);
colormap(blueColormap);
colorbar('TickLabelInterpreter', 'latex');
xlabel('Classe Prevista');
ylabel('Classe Verdadeira');
%title('Matriz de Confusão Após a Transferência de Aprendizado');
axis square;
grid off;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Incrustado', 'Limpo'});
set(gca, 'YTick', 1:2, 'YTickLabel', {'Incrustado', 'Limpo'});

% Anotar células com valores em porcentagem
textStrings = num2str(confMat_jda_percentage(:), '%.1f\\%%');
textStrings = strtrim(cellstr(textStrings));
[x, y] = meshgrid(1:2);
text(x(:), y(:), textStrings(:), 'HorizontalAlignment', 'center');

%% Função Auxiliar para Plotar Fronteira de Decis\~ao SVM
function svmBoundary(X, model)
    d = 0.01;
    [x1Grid, x2Grid] = meshgrid(min(X(:, 1)):d:max(X(:, 1)), ...
                                min(X(:, 2)):d:max(X(:, 2)));
    XGrid = [x1Grid(:), x2Grid(:)];
    [~, scores] = predict(model, XGrid);
    contour(x1Grid, x2Grid, reshape(scores(:, 2), size(x1Grid)), [0 0], 'k', 'LineWidth', 1.5);
end
