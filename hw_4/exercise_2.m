%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course:       Modern econometric and statistical learning
%               methods for quantitative asset management
%
% Instructor:   Prof. Dr. Marc Paolella, Urban Ulrych
%               University of Zurich
%
% Author:       Rino Beeli
%
% Date:         April 16th, 2020
% 
% Topic:        Homework 4 - Exercise 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all force; rng default;



% read file `m_sp500ret_3mtcm.txt`
sp500_3mtcm = readtable('data/m_sp500ret_3mtcm.txt', 'HeaderLines', 1, ...
    'Format', '%s	%s	%f	%s	%f', 'PreserveVariableNames', true);

% risk-free rate (monthly)
rf = sp500_3mtcm{:, '3mTCM'} / 100 / 12;

% S&P 500 returns
sp500 = sp500_3mtcm{:, 'sp500'};
sp500_ex = sp500 - rf;

clear sp500_3mtcm;



% read file `m_logret_10stocks.txt`
stocks = readtable('data/m_logret_10stocks.txt', 'Format', '%{MM/dd/yyyy}D %f %f %f %f %f %f %f %f %f %f');
stocks = stocks(~isnat(stocks{:, 1}), :); % remove empty rows
dates = stocks{:, 1}; % extract date
tickers = stocks.Properties.VariableNames(2:end);
stocks = stocks{:, 2:end}; % extract returns
[n, p] = size(stocks);
stocks_ex = stocks - rf;

% sanity check (same number of rows)
assert(size(sp500, 1) == size(stocks, 1));





% ----------------------------------------------
% a) Structured covariance matrix based on CAPM
% ----------------------------------------------

% structured covariance matrix single-factor model
F = CovMat1F(stocks_ex, sp500_ex);

% MLE of covariance matrix
S_MLE = cov(stocks_ex, 1); % MLE 1/n
S_MLE_diag = diag(S_MLE);


% plot covariance matrix, compare with sample covariance matrix
set(gcf, 'Position',[600, 400, 1400, 400])
tiledlayout(1, 3, 'Padding','compact')

nexttile; plotCovMat('MLE covariance matrix', S_MLE, 1:p, 1:p)
nexttile; plotCovMat('Structured cov. matrix (single-factor model)', F, 1:p, 1:p)




% ----------------------------------------------
% b) Shrink covariance matrix
% ----------------------------------------------

% expected returns
mu = mean(stocks_ex)';

% shrink sample covariance matrix
shrinkage = CovMatShrinkage(stocks_ex, F);
S_shrunk = shrinkage * F + (1 - shrinkage)*S_MLE;

% plot shrunk covariance matrix
nexttile; plotCovMat(sprintf('Shrunk covariance matrix (delta=%.3f)', shrinkage), S_shrunk, 1:p, 1:p)


% save to PDF
set(gcf, 'Units','inches')
pos = get(gcf, 'Position');   
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print(gcf, '-dpdf', '../report/2_cov_mat.pdf');





% compute efficient frontier portfolios using shrunk covariance matrix
Z = 50; % number of portfolios / ranks
ef_allowShorts = true;

ef_shrunk = EfficientFrontier(mu, S_shrunk, Z, ef_allowShorts);


% plot see c)





% ----------------------------------------------
% c) Efficient frontier using sample covariance
%    (no shrinkage)
% ----------------------------------------------

ef_MLE = EfficientFrontier(mu, S_MLE, Z, ef_allowShorts);


% plot efficient frontiers
figure
set(gcf, 'Position',  [900, 500, 700, 550])

% using sample covariance
plot(ef_MLE.Risk, ef_MLE.Return, 'r-', 'LineWidth',1); hold on;

% using 1-factor model shrinkage (S_shrunk space)
plot(ef_shrunk.Risk, ef_shrunk.Return, 'b-', 'LineWidth',1); hold on;

% % using 1-factor model shrinkage (S_MLE space)
% ef_shrunk_rets = zeros(Z, 1);
% ef_shrunk_volas = zeros(Z, 1);
% for i=1:Z
%     ef_shrunk_rets(i) = ef_shrunk.Weights(i, :) * mu;
%     ef_shrunk_volas(i) = sqrt(ef_shrunk.Weights(i, :) * S_MLE * ef_shrunk.Weights(i, :)');
% end
% plot(ef_shrunk_volas, ef_shrunk_rets, 'b--', 'LineWidth',1); hold on;

% assets
plot(sqrt(diag(S_MLE)), mu, 'k*', 'MarkerSize',4); hold on;
text(sqrt(diag(S_MLE))*1.02, mu, tickers, 'FontSize',10); hold on;

% legend({'EF MLE', 'EF Shrinkage S_{shrunk}', 'EF Shrinkage S_{MLE}', 'Assets'}, ...
%     'Location','southwest')

legend({'EF MLE', 'EF Shrinkage', 'Assets'}, 'Location','southwest')

xlabel('Standard deviation $\sigma$', 'Interpreter','latex')
ylabel('Expected excess return $\mu^e$', 'Interpreter','latex')
title('Mean-Variance Efficient Frontiers')
grid on
box on

    

% save to PDF
set(gcf, 'Units','inches')
pos = get(gcf, 'Position');   
set(gcf, 'PaperUnits','inches');
set(gcf, 'PaperSize', [pos(3) pos(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
print(gcf, '-dpdf', '../report/2_EFs.pdf');




% save variables used in exercise 3 to file `exercise_3.mat`
save('exercise_3.mat', 'stocks_ex', 'sp500_ex', 'tickers', 'mu', 'S_MLE', 'S_shrunk', 'F')




function plotCovMat(titleStr, covMat, xLabels, yLabels)
    % plots the covariance matrix covMat as heatmap.
    n = size(covMat, 1);
    imagesc(covMat) % plot the matrix
    set(gca, 'XTick', 1:n) % center x-axis ticks on bins
    set(gca, 'YTick', 1:n) % center y-axis ticks on bins
    set(gca, 'XTickLabelRotation', 90')
    set(gca, 'XTickLabel', xLabels)
    set(gca, 'YTickLabel', yLabels) 
    title(titleStr, 'FontSize', 14)
    colormap('jet')
    colorbar
end


