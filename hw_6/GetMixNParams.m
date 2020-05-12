function [wgts, means, alphas, betas, gammas] = GetMixNParams(x)
    % Extracts the MixN(3,2)-GARCH(1,1) model parameters
    % from the optimizer's parameters vector.
    % The last mean mu_k is chosen such that the
    % mean of the Gaussian mixture distribution is 0.
    % The weights of the G.M.D. sum to 1.
    
    % MixN weights
    wgts = x(1:3)';

    % MixN means
    means = zeros(3,1);
    means(1:end-1) = x(4:5);
    means(3) = -sum(wgts(1:end-1) .* means(1:end-1) ./ wgts(1:end-1)); % mean-zero
    
    % GARCH alphas (constant)
    alphas = x(6:8)';
    
    % GARCH betas (ARCH coef)
    betas = zeros(3,1);
    betas(1:2) = x(9:10);
    
    % GARCH gammas (GARCH coef)
    gammas = zeros(3,1);
    gammas(1:2) = x(11:end);
end