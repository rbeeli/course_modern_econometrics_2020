clc; clear; close all force; 



% ----------------------------------------------------------
% code based on example 1.11, page 37 of book Linear Models
% and Time-Series Analysis
% ----------------------------------------------------------

% set RNG seed (legacy implementation)
randn('state', 2);

cc = 5;  % T=40
%cc = 10; % T=80
T = 2*4*cc; % cc is cell count. So T is a multiple of 2*4
J = 2;
beta = [0 -5 3 5]';
%beta = [0 -1 3 4]';
dum1 = [ones(T/2, 1); zeros(T/2, 1)];
dum2 = 1 - dum1;
time = kron((1:4)', ones(cc, 1));
c3 = kron([1, 0]', time);
c4 = kron([0, 1]', time);
X = [dum1 dum2 c3 c4];
H = [1 -1 0 0;
     0 0 1 -1];
h = zeros(size(H, 1), 1);
XXinv = inv(X'*X);
HAHinv = inv(H*XXinv*H');




% -----------------------------------------------------------------------
% i. simulate regression many times and calculate power based on p-values
% -----------------------------------------------------------------------

sim = 200000;
sim_F_stats = zeros(sim, 1);
sim_F_pvals = zeros(sim, 1);
sim_LR = zeros(sim, 1);
sim_LR_stats = zeros(sim, 1);
sim_LR_pvals = zeros(sim, 1);
sim_betahat = zeros(sim, 4);

for i=1:sim
    % generate target variable with noise
    y = X*beta + 3*randn(T, 1);

    % unconstrained OLS regression
    betahat = XXinv*X'*y;
    sim_betahat(i, :) = betahat;
    yhat = X*betahat;
    res = y - yhat;
    Sbeta = sum(res.^2);
    sig2hat = Sbeta/(T - 4);

    % restriced OLS regression
    gammahat = betahat + XXinv*H'*HAHinv*(h - H*betahat);
    yhat2 = X*gammahat;
    res2 = y - yhat2;
    Sgamma = sum(res2.^2);
    sim_F_stats(i) = (Sgamma - Sbeta) / (J * sig2hat);

    % F-test for power computation
    sim_F_pvals(i) = 1 - fcdf(sim_F_stats(i), J, T-4);

    % likelihood ratio test for power computation
    sim_LR(i) = (Sgamma/Sbeta)^(-T/2);
    sim_LR_stats(i) = -2*log(sim_LR(i));
    sim_LR_pvals(i) = 1 - chi2cdf(sim_LR_stats(i), 2);

    % more precise F-test based on LR statistic
    % https://www.tandfonline.com/doi/abs/10.1080/03610926.2015.1096397
    %sim_LR_pvals(i) = 1 - fcdf((exp(sim_LR_stat(i)/T)-1)*(T-4)/J, J, T-4);
end


% results
fprintf('T=%i, beta=[%i, %i, %i, %i] \n', T, beta(1), beta(2), beta(3), beta(4))
fprintf('F-test  power at alpha=0.05: %.3f \n', mean(sim_F_pvals < 0.05))
fprintf('F-test  power at alpha=0.01: %.3f \n', mean(sim_F_pvals < 0.01))
fprintf('LR test power at alpha=0.05: %.3f (asymptotic Chi2(2)) \n', mean(sim_LR_pvals < 0.05))
fprintf('LR test power at alpha=0.01: %.3f (asymptotic Chi2(2)) \n\n', mean(sim_LR_pvals < 0.01))


% plots
figure
subplot(2,1,1)
histogram(sim_F_stats, 30, 'BinMethod','integers')
title(sprintf('F-test statistic distribution for T=%i, beta=[%i, %i, %i, %i]', T, beta(1), beta(2), beta(3), beta(4)))

subplot(2,1,2)
histogram(sim_LR_stats, 30, 'BinMethod','integers')
title(sprintf('LR-test statistic distribution for T=%i, beta=[%i, %i, %i, %i]', T, beta(1), beta(2), beta(3), beta(4)))



% ----------------------------------------------------
% ii. Bootstrap on OLS residuals for power calculation
% ----------------------------------------------------

% big outer for-loop
R = 2000;

% number of bootstrap replications
B = 500;

b_pvals_onesided = zeros(R, 1);
b_pvals_twosided = zeros(R, 1);

for r=1:R
    r_LR = NaN(B, 1);
    r_LR_stats = NaN(B, 1);
    r_LR_pvalues = NaN(B, 1);
    r_betas = NaN(B, 4);

    r_betahat = sim_betahat(r, :)';
    r_yhat = X*r_betahat;
    r_res = y - r_yhat;

    for b=1:B
        % sample OLS residuals with replacement
        b_u = r_res(randi(T, T, 1));

        % create synthetic Ys
        b_y = r_yhat + b_u;

        % unconstrained OLS regression
        b_betahat = XXinv*X'*b_y;
        b_yhat = X*b_betahat;
        b_res = b_y - b_yhat;
        b_Sbeta = sum(b_res.^2);
        r_betas(b, :) = b_betahat;

        % restriced OLS regression
        b_gammahat = b_betahat + XXinv*H'*HAHinv*(h - H*b_betahat);
        b_yhat2 = X*b_gammahat;
        b_res2 = b_y - b_yhat2;
        b_Sgamma = sum(b_res2.^2);

        % likelihood ratio
        r_LR(b) = (b_Sgamma/b_Sbeta)^(-T/2);
        r_LR_stats(b) = -2*log(r_LR(b));
    end

    % two-sided test: beta1 != beta2 and beta3 != beta4
    pvalue_test = 1 - r_LR < sim_LR(r);
    b_pvals_twosided(r) = mean(pvalue_test);

    % one-sided test: beta1 > beta2 and beta3 > beta4
    onsided_test = (r_betas(:, 1) > r_betas(:, 2)) & (r_betas(:, 3) < r_betas(:, 4));
    b_pvals_onesided(r) = mean(pvalue_test & onsided_test);
end



two_sid_pwr05 = mean(b_pvals_twosided < 0.05);
two_sid_pwr01 = mean(b_pvals_twosided < 0.01);

one_sid_pwr05 = mean(b_pvals_onesided < 0.05);
one_sid_pwr01 = mean(b_pvals_onesided < 0.01);

fprintf('Two-sided bootstrap LR power at alpha=0.05: %.3f \n', two_sid_pwr05)
fprintf('Two-sided bootstrap LR power at alpha=0.01: %.3f \n\n', two_sid_pwr01)

fprintf('One-sided bootstrap LR power at alpha=0.05: %.3f \n', one_sid_pwr05)
fprintf('One-sided bootstrap LR power at alpha=0.01: %.3f \n\n', one_sid_pwr01)

% how much does the power improve when moving from two- to one-sided tests?
fprintf('One-sided power improvement at alpha=0.05: %.2f%% \n', 100*(one_sid_pwr05/two_sid_pwr05-1))
fprintf('One-sided power improvement at alpha=0.01: %.2f%% \n', 100*(one_sid_pwr01/two_sid_pwr01-1))
