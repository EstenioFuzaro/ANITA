%% ------------------------------------------------------------------------- %%
% Transfer Learning - Classification with Thresholds for HE1 and HE2
% --------------------------------------------------------------------------

clc;
clear;
close all;

%---------------------------------------------------------------------------

% Configurações de Estilo
nfonte = 20;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesXMinorGrid', 'on', 'defaultAxesXMinorGridMode', 'manual');
set(groot, 'defaultAxesYMinorGrid', 'on', 'defaultAxesYMinorGridMode', 'manual');
set(groot, 'defaultAxesZMinorGrid', 'on', 'defaultAxesZMinorGridMode', 'manual');

%---------------------------------------------------------------------------

% Definir Parâmetros Estocásticos
rng_stream = RandStream('mt19937ar', 'Seed', 30081984);
RandStream.setGlobalStream(rng_stream);

%---------------------------------------------------------------------------

% Carregar Features
load('features_TC1.mat'); % Features do HE1
X1_source = X_1(1:199);
X2_source = X_2(1:199);

load('features_TC2.mat'); % Features do HE2
X1_target = X_1(1:39);
X2_target = X_2(1:39);

% Combinar Dados do HE1 e HE2
X_source = [X1_source(:), X2_source(:)];
X_target = [X1_target(:), X2_target(:)];

% Rótulos para HE1 (HE1)
Y_source = [ones(99, 1); -ones(100, 1)];

% Rótulos para HE2 (HE2)
Y_target = [ones(19, 1); -ones(20, 1)];

%---------------------------------------------------------------------------
% Parte 1: Sem Normalização
%---------------------------------------------------------------------------

% Treinamento
X_train_clean = X_source(1:99, :);
Y_train_clean = Y_source(1:99);
X_train_fouled = X_source(100:end, :);
Y_train_fouled = Y_source(100:end);
X_train = [X_train_clean; X_train_fouled];
Y_train = [Y_train_clean; Y_train_fouled];

% Validação
X_val = [X_source; X_target];
Y_val = [Y_source; Y_target];

% Treinar SVM
SVMModel = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'Standardize', false);
Y_pred_val = predict(SVMModel, X_val);

% Visualizar Fronteira de Decisão
visualizeDecisionBoundary(X_source, X_target, X_val, Y_pred_val, 'HE1 + HE2 (Original)', nfonte, SVMModel);
 % Configurações do Gráfico
    xlabel('$\mathcal{X}_1$', 'Interpreter', 'latex', 'FontSize', nfonte);
    ylabel('$\mathcal{X}_2$', 'Interpreter', 'latex', 'FontSize', nfonte);
    legend('Location', 'Best', 'Interpreter', 'latex');
    grid on;
    hold off;

% Matriz de Confusão
confMat_val = confusionmat(Y_val, Y_pred_val);
plotConfMatrix(confMat_val, 'Confusion Matrix - Validation (Original)');

% Visualização com `gplotmatrix` (Sem Normalização)
% Combine features and labels
features_combined = [X_source; X_target];
tc_labels = [ones(size(X_source, 1), 1); 2 * ones(size(X_target, 1), 1)];
group_tc = categorical(tc_labels, [1 2], {'HE1 (Source)', 'HE2 (Target)'});
condition_labels = [repmat({'Clean'}, 99, 1); repmat({'Fouled'}, 100, 1);
                    repmat({'Clean'}, 19, 1); repmat({'Fouled'}, 20, 1)];
group_condition = categorical(condition_labels, {'Clean', 'Fouled'});
combined_labels = strcat(cellstr(group_tc), '-', cellstr(group_condition));
xnames = {'$\mathcal{X}_1$', '$\mathcal{X}_2$'};

% Define colors and marker symbols
colors = [0 0 1; 1 0 0; 0 1 0; 1 0 1]; % Example RGB colors for groups
markerSymbols = 'osd^'; % Circle, square, diamond, triangle

% Full-screen figure
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Create scatter plot matrix
[h, ax] = gplotmatrix(features_combined, [], combined_labels, colors, markerSymbols, [], [], 'grpbars', xnames, xnames);

% Customize markers (black edges)
for i = 1:size(h, 1)
    for j = 1:size(h, 2)
        for k = 1:size(h, 3) % Iterate over groups
            % Check if the handle is a scatter plot
            if ~isempty(h(i, j, k)) && isgraphics(h(i, j, k), 'line')
                h(i, j, k).MarkerFaceColor = h(i, j, k).Color; % Match fill color to marker color
                h(i, j, k).MarkerEdgeColor = 'k'; % Set edge color to black
            end
        end
    end
end

% Add title and adjust font sizes
%title('Feature Distribution (Customized)', 'Interpreter', 'latex');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

%---------------------------------------------------------------------------
% Normalização
%---------------------------------------------------------------------------

% Normalizar Dados
X1_source = zscore(X1_source);
X2_source = zscore(X2_source);
X1_target = zscore(X1_target);
X2_target = zscore(X2_target);

X_source = [X1_source(:), X2_source(:)];
X_target = [X1_target(:), X2_target(:)];

% Treinamento com Dados Normalizados
X_train_clean = X_source(1:99, :);
X_train_fouled = X_source(100:end, :);
X_train = [X_train_clean; X_train_fouled];

% Validação com Dados Normalizados
X_val = [X_source; X_target];

% Treinar SVM com Dados Normalizados
SVMModel_norm = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'Standardize', false);
Y_pred_val_norm = predict(SVMModel_norm, X_val);

% Visualizar Fronteira de Decisão (Normalizado)
visualizeDecisionBoundary(X_source, X_target, X_val, Y_pred_val_norm, ...
                          'HE1 + HE2 (Normalized)', nfonte, SVMModel_norm);
 % Configurações do Gráfico
    xlabel('$\mathcal{X}_1$', 'Interpreter', 'latex', 'FontSize', nfonte);
    ylabel('$\mathcal{X}_2$', 'Interpreter', 'latex', 'FontSize', nfonte);
    legend('Location', 'Best', 'Interpreter', 'latex');
    grid on;
    hold off;

% Matriz de Confusão (Normalizada)
confMat_val_norm = confusionmat(Y_val, Y_pred_val_norm);
plotConfMatrix(confMat_val_norm, 'Confusion Matrix - Validation (Normalized)');

% Combine normalized features and labels
features_combined_norm = [X_source; X_target];

% Full-screen figure
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

% Define colors and marker symbols
colors = [0 0 1; 1 0 0; 0 1 0; 1 0 1]; % Example RGB colors for groups
markerSymbols = 'osd^'; % Circle, square, diamond, triangle

% Create scatter plot matrix
[h, ax] = gplotmatrix(features_combined_norm, [], combined_labels, colors, markerSymbols, [], [], 'grpbars', xnames, xnames);

% Customize markers (black edges)
for i = 1:size(h, 1)
    for j = 1:size(h, 2)
        for k = 1:size(h, 3) % Iterate over groups
            % Check if the handle is a scatter plot
            if ~isempty(h(i, j, k)) && isgraphics(h(i, j, k), 'line')
                h(i, j, k).MarkerFaceColor = h(i, j, k).Color; % Match fill color to marker color
                h(i, j, k).MarkerEdgeColor = 'k'; % Set edge color to black
            end
        end
    end
end

% Add title and adjust font sizes
%title('Feature Distribution (Normalized)', 'Interpreter', 'latex');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

% ------------------------------------------------------------------------- 
% Transfer Learning using JDA
% -------------------------------------------------------------------------

% Source data (HE1)
Xs = [X1_source(:), X2_source(:)]; % Combine X1 and X2 for HE1
Ys = Y_source;                     % Labels for HE1

% Target data (HE2)
Xt = [X1_target(:), X2_target(:)]; % Combine X1 and X2 for HE2
Yt = Y_target;                     % Labels for HE2 (for evaluation only)

% Parameters for JDA
lambda = 1e-3;          % Regularization parameter
dim = 2;                % Number of dimensions after projection
kernel_type = 'linear'; % Kernel type ('linear' or 'rbf')

% -------------------------------------------------------------------------
% Perform JDA Transformation
% -------------------------------------------------------------------------
[Zs, Zt] = JDA(Xs, Xt, Ys, lambda, dim, kernel_type);

% Combine Transformed Data for Visualization
Z = [Zs; Zt];
Y = [Ys; Yt];

% -------------------------------------------------------------------------
% Scatter Plot with Decision Boundary (JDA)
% -------------------------------------------------------------------------
% Train SVM on Transformed Source Data
SVMModel_JDA = fitcsvm(Zs, Ys, 'KernelFunction', 'linear', 'Standardize', false);

% Predict on Transformed Target Data
Y_pred_JDA = predict(SVMModel_JDA, Z);

% Scatter Plot with Decision Boundary
visualizeDecisionBoundary(Zs, Zt, Z, Y_pred_JDA, 'JDA Transformed Features', nfonte, SVMModel_JDA);
 % Configurações do Gráfico
    xlabel('$\mathcal{Z}_1$', 'Interpreter', 'latex', 'FontSize', nfonte);
    ylabel('$\mathcal{Z}_2$', 'Interpreter', 'latex', 'FontSize', nfonte);
    legend('Location', 'Best', 'Interpreter', 'latex');
    grid on;
    hold off;

% -------------------------------------------------------------------------
% Confusion Matrix for JDA
% -------------------------------------------------------------------------
confMat_JDA = confusionmat(Y, Y_pred_JDA);
plotConfMatrix(confMat_JDA, 'Confusion Matrix - JDA (HE2)');

% -------------------------------------------------------------------------
% Visualize Transformed Features with `gplotmatrix`
% -------------------------------------------------------------------------

% Create Labels for Visualization
tc_labels = [ones(size(Zs, 1), 1); 2 * ones(size(Zt, 1), 1)];
group_tc = categorical(tc_labels, [1, 2], {'HE1 (Source)', 'HE2 (Target)'});
condition_labels = [repmat({'Clean'}, 99, 1); repmat({'Fouled'}, 100, 1);
                    repmat({'Clean'}, 19, 1); repmat({'Fouled'}, 20, 1)];
group_condition = categorical(condition_labels, {'Clean', 'Fouled'});
combined_labels = strcat(cellstr(group_tc), '-', cellstr(group_condition));
xnames = {'$\mathcal{Z}_1$', '$\mathcal{Z}_2$'};

% Define colors and marker symbols
colors = [0 0 1; 1 0 0; 0 1 0; 1 0 1]; % Example RGB colors for groups
markerSymbols = 'osd^'; % Circle, square, diamond, triangle

% Plot the Transformed Data
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % Fullscreen
[h, ax] = gplotmatrix(Z, [], combined_labels, colors, markerSymbols, [], [], 'grpbars', xnames, xnames);

% Customize markers (black edges)
for i = 1:size(h, 1)
    for j = 1:size(h, 2)
        for k = 1:size(h, 3) % Iterate over groups
            % Check if the handle is a scatter plot
            if ~isempty(h(i, j, k)) && isgraphics(h(i, j, k), 'line')
                h(i, j, k).MarkerFaceColor = h(i, j, k).Color; % Match fill color to marker color
                h(i, j, k).MarkerEdgeColor = 'k'; % Set edge color to black
            end
        end
    end
end

% Add title and adjust font sizes
%title('Feature Distribution after JDA', 'Interpreter', 'latex');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

%---------------------------------------------------------------------------
% Funções Auxiliares
%------------------------------------------------------------------- 
function plotConfMatrix(confMat, titleText)
    figure;
    numColors = 100;
    blueColormap = [linspace(0.678, 0, numColors)', ...
                    linspace(0.847, 0, numColors)', ...
                    linspace(0.902, 0.545, numColors)'];
    confMat_percentage = (confMat ./ sum(confMat, 2)) * 100; % Porcentagens por classe
    confMat_percentage(isnan(confMat_percentage)) = 0; % Corrigir possíveis NaN
    imagesc(confMat_percentage);
    colormap(blueColormap);
    colorbar('TickLabelInterpreter', 'latex');
    xlabel('Predicted Class', 'Interpreter', 'latex');
    ylabel('True Class', 'Interpreter', 'latex');
    %title(titleText, 'Interpreter', 'latex');
    axis square;
    grid off;
    set(gca, 'XTick', 1:2, 'XTickLabel', {'Fouled', 'Clean'}, 'TickLabelInterpreter', 'latex');
    set(gca, 'YTick', 1:2, 'YTickLabel', {'Fouled', 'Clean'}, 'TickLabelInterpreter', 'latex');

    % Adiciona texto com porcentagens
    textStrings = num2str(confMat_percentage(:), '%.1f\\%%');
    textStrings = strtrim(cellstr(textStrings));
    [x, y] = meshgrid(1:2);
    for i = 1:numel(textStrings)
        textColor = 'black';
        if confMat(i) > 0 && x(i) == y(i) % Destaque na diagonal
            textColor = 'white';
        end
        text(x(i), y(i), textStrings{i}, 'HorizontalAlignment', 'center', ...
            'Color', textColor, 'FontWeight', 'bold', 'FontSize', 14);
    end
end

%------------------------------------------------------------------- 
function visualizeDecisionBoundary(X_source, X_target, X_val, Y_pred_val, titleText, nfonte, model)
    % Function for visualizing the decision boundary
    % and points from source (HE1) and target (HE2)
    
    figure;
    hold on;

    % Plot Data from Source (HE1) with filled markers and black edges
    scatter(X_source(1:99, 1), X_source(1:99, 2), 80, 'b', 'o', ...
        'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'HE1 (Source)-Clean');
    scatter(X_source(100:end, 1), X_source(100:end, 2), 80, 'r', '^', ...
        'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'HE1 (Source)-Fouled');

    % Plot Data from Target (HE2) with filled markers and black edges
    scatter(X_target(1:19, 1), X_target(1:19, 2), 80, 'g', 's', ...
        'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'HE2 (Target)-Clean');
    scatter(X_target(20:end, 1), X_target(20:end, 2), 80, 'm', 'd', ...
        'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'HE2 (Target)-Fouled');

    % Uncomment to Plot Detected Fouled Points
    % fouled_detected_idx = find(Y_pred_val == -1); % Points classified as fouled
    % scatter(X_val(fouled_detected_idx, 1), X_val(fouled_detected_idx, 2), ...
    %         100, 'k', '+', 'LineWidth', 2, 'DisplayName', 'Detected Fouled');

    % Plot Decision Boundary
    svmBoundary([X_source; X_target], model);

    % Add title and adjust font size
    %title(titleText, 'FontSize', nfonte);
    legend('Location', 'bestoutside');
    hold off;
end



%------------------------------------------------------------------- 
function svmBoundary(X, model)
    d = 0.01;
    [x1Grid, x2Grid] = meshgrid(min(X(:, 1)):d:max(X(:, 1)), ...
                                min(X(:, 2)):d:max(X(:, 2)));
    XGrid = [x1Grid(:), x2Grid(:)];
    [~, scores] = predict(model, XGrid);
    contour(x1Grid, x2Grid, reshape(scores(:, 2), size(x1Grid)), [0 0], 'k', 'LineWidth', 1.5, ...
            'DisplayName', 'Decision Boundary');
end
