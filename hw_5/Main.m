%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course:       Modern econometric and statistical learning
%               methods for quantitative asset management
%
% Instructor:   Prof. Dr. Marc Paolella, Urban Ulrych
%               University of Zurich
%
% Author:       Rino Beeli
%
% Date:         April 28th, 2020
% 
% Topic:        Homework 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all force; rng default;
addpath('lib')



% AR(4) model coefficients
a_n = 4; a = [repmat([0.4 -0.3 0.2], a_n, 1) linspace(-0.1, -0.6, a_n)'];

% MA(3) model coefficients
b_n = 3; b = [repmat(-0.1368, b_n, 1) linspace(-0.8673, -0.1956, b_n)' repmat(0.0046, b_n, 1)];


T = 30;                                   % sample size
rep_spa = 250;                            % number of simulations using SPA
rep_lasso = 2000;                         % number of simulations using LASSO
p_max = floor(1.5*sqrt(T));               % max. p* - used for c in ButPao(Y, X, c)
B = 150;                                  % number of SPA bootstrap samples
X = [ones(T,1), (1:T)'];                  % constant + time trend design matrix
alpha_c = linspace(0.175, 0.025, p_max);  % significance levels c used in ButPao(Y,X,c)
alpha_coefs = 0.05;                       % significance level(s) LR test (scalar or vector, bonus question!)
%alpha_coefs = linspace(0.175, 0.025, p_max);




% simulate every model and save plots to file
models = cell(0, 2); % cell array of {model type, coefficients}
for i=1:a_n, models(end+1, :) = {'AR', a(i, :)}; end
for i=1:b_n, models(end+1, :) = {'MA', b(i, :)}; end
for i=1:(a_n + b_n)
    % get model type (MA, AR) and its coefficients
    [type, m_coefs] = models{i, :};
    
    % DGP function - simulates AR/MA time series
    if strcmp(type, 'AR')
        dgp = @(T) armasim(T, 1, m_coefs, []); % AR
    else
        dgp = @(T) armasim(T, 1, [], m_coefs); % MA
    end
    
    % simulate AR/MA data for LASSO and SPA
    Y_lasso = zeros(rep_lasso, T);
    Y_spa = zeros(rep_spa, T);
    for j=1:rep_lasso, Y_lasso(j, :) = dgp(T); end
    for j=1:rep_spa, Y_spa(j, :) = dgp(T); end

    
    % <LASSO>
    tic; [lasso_num_coefs, lasso_coefs, lasso_coefs_sig] = AR_LASSO(X, Y_lasso, p_max);
    
    % unique filename / identifier
    model_desc = sprintf('%s[%s] T=%i rep=%i pmax=%i', type, join(string(m_coefs), ','), T, rep_lasso, p_max);
    filename = sprintf('LASSO %s', model_desc);
    disp(model_desc)
    
    % plot and save PDFs
    close all force; GenPlots('LASSO', filename, model_desc, lasso_num_coefs, lasso_coefs, lasso_coefs_sig, p_max);
    fclose(fopen(sprintf('plots/%s_runtime=%.1f.txt', filename, toc), 'wb'));
    % </LASSO>
    
    
    % <SPA>
    tic; [spa_num_coefs, spa_coefs, spa_coefs_sig, LR, LRboot] = AR_SPA(X, Y_spa, B, alpha_c, alpha_coefs);
    
    % unique filename / identifier
    model_desc = sprintf('%s[%s] T=%i rep=%i B=%i pmax=%i', type, join(string(m_coefs), ','), T, rep_spa, B, p_max);
    filename = sprintf('SPA %s', model_desc);
    disp(model_desc)
    
    % plot and save PDFs
    close all force; GenPlots('SPA', filename, model_desc, spa_num_coefs, spa_coefs, spa_coefs_sig, p_max, LRboot);
    fclose(fopen(sprintf('plots/%s_runtime=%.1f.txt', filename, toc), 'wb'));
    % </SPA>
end








