function F = CovMat1F(assetRets, mktRets)
    % Computes the structured covariance matrix based on the
    % single-factor model CAPM. Regresses all asset returns
    % onto the supplied market returns.
    % It is possible to vectorize this implementation for
    % better performance.
    
    [N, K] = size(assetRets);
    D = zeros(K, K);
    betas = NaN(K, 2);
    
    for i=1:K
        % regress assets onto market
        y = assetRets(:, i);
        X = [ones(N, 1) mktRets]; % add intercept (ones)

        % perform regression
        [beta, ~, res] = regress(y, X);
        betas(i, :) = beta;

        % residuals variance
        D(i, i) = var(res, 1);
    end

    % structured covariance matrix
    F = var(mktRets, 1) * (betas(:, 2)*betas(:, 2)') + D;
end
