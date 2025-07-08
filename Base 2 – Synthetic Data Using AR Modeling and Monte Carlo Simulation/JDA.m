function [Zs, Zt] = JDA(Xs, Xt, Ys, lambda, dim, kernel_type)
    % Inputs:
    % Xs: Source domain data (m x d)
    % Xt: Target domain data (n x d)
    % Ys: Source domain labels (m x 1)
    % lambda: Regularization parameter (scalar)
    % dim: Number of dimensions after mapping (scalar)
    % kernel_type: Type of kernel ('linear', 'rbf', etc.)

    % Combine source and target data
    X = [Xs; Xt];
    n = size(Xs, 1);
    m = size(Xt, 1);
    C = length(unique(Ys));  % Number of classes
    X = X';

    % Compute Kernel matrix
    if strcmp(kernel_type, 'linear')
        K = X' * X;
    elseif strcmp(kernel_type, 'rbf')
        sigma = mean(mean(pdist2(X', X', 'euclidean').^2));
        K = exp(-pdist2(X', X', 'euclidean').^2 / (2 * sigma));
    else
        error('Unknown kernel_type');
    end

    % Compute MMD matrix
    e = [1/n * ones(n, 1); -1/m * ones(m, 1)];
    M = e * e' * C;
    for c = unique(Ys)'
        e = zeros(n + m, 1);
        e(Ys == c) = 1/length(find(Ys == c));
        M = M + e * e';
    end
    M = M / norm(M, 'fro');

    % Centering matrix
    H = eye(n + m) - 1/(n + m) * ones(n + m, n + m);

    % Compute JDA
    A = (K * M * K' + lambda * eye(n + m)) \ (K * H * K');
    [V, ~] = eigs(A, dim);

    % Project source and target data
    Z = V' * K;
    Zs = Z(:, 1:n)';
    Zt = Z(:, n+1:end)';
end
