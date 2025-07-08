%% ------------------------------------------------------------------------- %%
% Transfer Learning - Classification with Thresholds for Shell-Tube and Plates
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
rng_stream = RandStream('mt19937ar', 'Seed', 16021971);
RandStream.setGlobalStream(rng_stream);

%---------------------------------------------------------------------------

n_exp = 5;
n_sim = 45;
n_fou = 50;
n_test = n_exp+n_sim+n_fou-1;
n_clean_test = n_exp+n_sim-1;

% Carregar Features
load('features_CascoTubo.mat'); % Features do Shell-Tube
X1_source = X_1(1:n_test);
X2_source = X_2(1:n_test);

load('features_PlacaPlana.mat'); % Features do Plates
X1_target = X_1(1:n_test);
X2_target = X_2(1:n_test);

% Combinar Dados do Shell-Tube e Plates
X_source = [X1_source(:), X2_source(:)];
X_target = [X1_target(:), X2_target(:)];

% Rótulos para Shell-Tube (Source)
Y_source = [ones(n_clean_test, 1); -ones(n_fou, 1)];

% Rótulos para Plates (Target)
Y_target = [ones(n_clean_test, 1); -ones(n_fou, 1)];

%---------------------------------------------------------------------------

%%  Parte 1: Treinamento de SVM sem normalização
%---------------------------------------------------------------------------

% Conjuntos de Treinamento e Validação
X_train = X_source;
Y_train = Y_source;

X_val = [X_source; X_target];
Y_val = [Y_source; Y_target];

% Treinar SVM
SVMModel = fitcsvm(X_train, Y_train, 'KernelFunction', 'linear', 'Standardize', false);
Y_pred_val = predict(SVMModel, X_val);

%  Visualizar Fronteira de Decisão
visualizeDecisionBoundary(X_source, X_target, X_val, Y_pred_val, 'Shell-Tube + Plates (Original)', nfonte, SVMModel);
% Add axis labels with LaTeX formatting
xlabel('$\mathcal{X}_1$', 'Interpreter', 'latex', 'FontSize', nfonte);
ylabel('$\mathcal{X}_2$', 'Interpreter', 'latex', 'FontSize', nfonte);

%  Matriz de Confusão
confMat_val = confusionmat(Y_val, Y_pred_val);
plotConfMatrix(confMat_val, 'Confusion Matrix - Validation (Original)');

% Visualização com `gplotmatrix` (Sem Normalização)
% Combine features and labels
features_combined = [X_source; X_target];
tc_labels = [ones(size(X_source, 1), 1); 2 * ones(size(X_target, 1), 1)];
group_tc = categorical(tc_labels, [1 2], {'Shell-Tube (Source)', 'Plates (Target)'});
condition_labels = [repmat({'Clean'}, n_clean_test, 1); repmat({'Fouled'}, n_fou, 1);
    repmat({'Clean'}, n_clean_test, 1); repmat({'Fouled'}, n_fou, 1)];
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

%%  Parte 2: Normalização dos Dados
%---------------------------------------------------------------------------

% Normalizar Dados
X1_source = zscore(X1_source);
X2_source = zscore(X2_source);
X1_target = zscore(X1_target);
X2_target = zscore(X2_target);

X_source = [X1_source(:), X2_source(:)];
X_target = [X1_target(:), X2_target(:)];

% Treinar SVM com Dados Normalizados
SVMModel_norm = fitcsvm(X_source, Y_source, 'KernelFunction', 'linear', 'Standardize', false);
Y_pred_val_norm = predict(SVMModel_norm, [X_source; X_target]);

%  Visualizar Fronteira de Decisão (Normalizado)
visualizeDecisionBoundary(X_source, X_target, [X_source; X_target], Y_pred_val_norm, ...
    'Shell-Tube + Plates (Normalized)', nfonte, SVMModel_norm);
% Add axis labels with LaTeX formatting
xlabel('$\mathcal{X}_1$', 'Interpreter', 'latex', 'FontSize', nfonte);
ylabel('$\mathcal{X}_2$', 'Interpreter', 'latex', 'FontSize', nfonte);

%  Matriz de Confusão (Normalizada)
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
%---------------------------------------------------------------------------

%%  Parte 3: Transfer Learning com Joint Distribution Adaptation (JDA)
%---------------------------------------------------------------------------

% Fonte (Shell-Tube) e Alvo (Plates)
Xs = X_source;  % Source data
Ys = Y_source;  % Source labels
Xt = X_target;  % Target data
Yt = Y_target;  % Target labels

% Parâmetros do JDA
lambda = 1e-3;
dim = 2;
kernel_type = 'linear';

%  Aplicar JDA
[Zs, Zt] = JDA(Xs, Xt, Ys, lambda, dim, kernel_type);
Z = [Zs; Zt];
Y = [Ys; Yt];

%  Treinar SVM no Domínio Adaptado
SVMModel_JDA = fitcsvm(Zs, Ys, 'KernelFunction', 'linear', 'Standardize', false);
Y_pred_JDA = predict(SVMModel_JDA, Z);

%  Visualizar Fronteira de Decisão com JDA
visualizeDecisionBoundary(Zs, Zt, Z, Y_pred_JDA, 'JDA Transformed Features', nfonte, SVMModel_JDA);
% Add axis labels with LaTeX formatting
xlabel('$\mathcal{Z}_1$', 'Interpreter', 'latex', 'FontSize', nfonte);
ylabel('$\mathcal{Z}_2$', 'Interpreter', 'latex', 'FontSize', nfonte);

%  Matriz de Confusão para JDA
confMat_JDA = confusionmat(Y, Y_pred_JDA);
plotConfMatrix(confMat_JDA, 'Confusion Matrix - JDA (Plates)');

% Create Labels for Visualization
tc_labels = [ones(size(Zs, 1), 1); 2 * ones(size(Zt, 1), 1)];
group_tc = categorical(tc_labels, [1, 2], {'Shell-Tube (Source)', 'Plates (Target)'});
condition_labels = [repmat({'Clean'}, n_clean_test, 1); repmat({'Fouled'}, n_fou, 1);
    repmat({'Clean'}, n_clean_test, 1); repmat({'Fouled'}, n_fou, 1)];
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
% and points from source (Shell-Tube) and target (Plates)
n_exp = 5;
n_sim = 45;
n_fou = 50;
n_test = n_exp+n_sim+n_fou-1;
n_clean_test = n_exp+n_sim-1;
figure;
box on;
hold on;

% Plot Data from Source (Shell-Tube) with filled markers and black edges
scatter(X_source(1:n_clean_test, 1), X_source(1:n_clean_test, 2), 80, 'b', 'o', ...
    'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Shell-Tube (Source)-Clean');
scatter(X_source(n_fou:end, 1), X_source(n_fou:end, 2), 80, 'r', '^', ...
    'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Shell-Tube (Source)-Fouled');

% Plot Data from Target (Plates) with filled markers and black edges
scatter(X_target(1:n_clean_test, 1), X_target(1:n_clean_test, 2), 80, 'g', 's', ...
    'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Plates (Target)-Clean');
scatter(X_target(n_fou:end, 1), X_target(n_fou:end, 2), 80, 'm', 'd', ...
    'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Plates (Target)-Fouled');

% Uncomment to Plot Detected Fouled Points
% fouled_detected_idx = find(Y_pred_val == -1); % Points classified as fouled
% scatter(X_val(fouled_detected_idx, 1), X_val(fouled_detected_idx, 2), ...
%         100, 'k', '+', 'LineWidth', 2, 'DisplayName', 'Detected Fouled');

% Plot Decision Boundary
svmBoundary([X_source; X_target], model);

% Add title and adjust font size
%title(titleText, 'FontSize', nfonte);
legend('Location', 'best');
hold off;
end

function svmBoundary(X, model)
n_exp = 5;
n_sim = 45;
n_fou = 50;
n_test = n_exp+n_sim+n_fou-1;
n_clean_test = n_exp+n_sim-1;
% Define the number of grid points (limits mesh size)
num_grid_points = 100;  % Set to a reasonable value to avoid memory overload

% Compute reasonable bounds while ignoring outliers
x1_min = prctile(X(:,1), 1); % Avoid extreme outliers (1st percentile)
x1_max = prctile(X(:,1), n_test);

x2_min = prctile(X(:,2), 1);
x2_max = prctile(X(:,2), n_test);

% Create grid with fixed number of points
x1_range = linspace(x1_min, x1_max, num_grid_points);
x2_range = linspace(x2_min, x2_max, num_grid_points);
[x1Grid, x2Grid] = meshgrid(x1_range, x2_range);

% Prepare grid for SVM prediction
XGrid = [x1Grid(:), x2Grid(:)];
[~, scores] = predict(model, XGrid);

% Plot decision boundary
contour(x1Grid, x2Grid, reshape(scores(:,2), size(x1Grid)), [0 0], 'k', 'LineWidth', 1.5, 'DisplayName', 'Decision Boundary');
end

