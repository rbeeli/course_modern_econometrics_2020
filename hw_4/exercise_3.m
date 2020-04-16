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
% Topic:        Homework 4 - Exercise 3.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all force; rng default;


% load data as calculated in exercise 2
load('exercise_3.mat')

% number of ranks
Z = 50;



% -------------------------
% Mean-Variance portfolios
% -------------------------

% allow short-selling?
allowShorts = false;

% Mean-Variance efficient frontier based on MLE
ef_MV_MLE = EfficientFrontier(mu, S_MLE, Z, allowShorts);

% Mean-Variance efficient frontier based on Shrinkage
ef_MV_shrink = EfficientFrontier(mu, S_shrunk, Z, allowShorts);



% -------------------------
% Portfolio resampling
% -------------------------

[T, K] = size(stocks_ex);

% N = number of draws/simulations
for N=[20 100]

    % store average resampled weight per rank
    avg_wgts_RE = zeros(Z, K);
    avg_wgts_RE_shrink = zeros(Z, K);

    for i=1:N
        % parametric bootstrapping of returns

        % draw T samples
        r_hat = mvnrnd(mu, S_MLE, T);
        r_hat_shrink = mvnrnd(mu, S_shrunk, T); % shrinkage

        % expected returns of samples
        mu_hat = mean(r_hat)';
        mu_hat_shrink = mean(r_hat_shrink)';

        % MLE covariance matrix
        S_r_hat = cov(r_hat, 1);
        S_r_hat_shrink = cov(r_hat_shrink, 1);

        % compute efficient frontier portfolios along Z ranks
        % starts at the GMV portfolio until maximum expected return
        efs_RE(i) = EfficientFrontier(mu_hat, S_r_hat, Z, allowShorts);
        efs_RE_shrink(i) = EfficientFrontier(mu_hat_shrink, S_r_hat_shrink, Z, allowShorts);

        % store average resampled weight per rank
        avg_wgts_RE = avg_wgts_RE + 1/N*efs_RE(i).Weights;
        avg_wgts_RE_shrink = avg_wgts_RE_shrink + 1/N*efs_RE_shrink(i).Weights;

        fprintf('Simulation i=%i of %i\n', i, N)
    end

    % resampled portfolios returns
    rets_RE = avg_wgts_RE * mu;
    rets_RE_shrink = avg_wgts_RE_shrink * mu;

    % resampled portfolios volatilities
    vola_RE = zeros(Z, 1);
    vola_RE_shrink = zeros(Z, 1);
    for i=1:Z
        vola_RE(i) = sqrt(avg_wgts_RE(i, :) * S_MLE * avg_wgts_RE(i, :)');
        vola_RE_shrink(i) = sqrt(avg_wgts_RE_shrink(i, :) * S_shrunk * avg_wgts_RE_shrink(i, :)');
    end




    % plot efficient frontiers
    figure
    set(gcf, 'Position',  [800, 200, 1000, 880])
    tiledlayout(1, 1, 'Padding','compact', 'TileSpacing','Compact')
    nexttile
    
    % all resampled portfolios using MLE covariance matrix
    for i=1:N
        for j=1:Z
            vola_j = sqrt(efs_RE(i).Weights(j, :) * S_MLE * efs_RE(i).Weights(j, :)');
            expRet_j = efs_RE(i).Weights(j, :) * mu;
            plt = plot(vola_j, expRet_j, 'm*', 'MarkerSize',0.5, 'HandleVisibility','off'); hold on;
        end
    end

    % all resampled portfolios using shrunk covariance matrix
    for i=1:N
        for j=1:Z
            vola_j = sqrt(efs_RE_shrink(i).Weights(j, :) * S_shrunk * efs_RE_shrink(i).Weights(j, :)');
            expRet_j = efs_RE_shrink(i).Weights(j, :) * mu;
            plt = plot(vola_j, expRet_j, 'c*', 'MarkerSize',0.5, 'HandleVisibility','off'); hold on;
        end
    end
    
    % MV MLE efficient frontier
    hold on; plot(ef_MV_MLE.Risk, ef_MV_MLE.Return, 'r-', 'LineWidth',1.5)

    % MV shrinkage efficient frontier
    hold on; plot(ef_MV_shrink.Risk, ef_MV_shrink.Return, 'b-', 'LineWidth',1.5)

    % averaged resampled efficient frontier
    hold on; plot(vola_RE, rets_RE, 'm-', 'LineWidth',1.5)

    % averaged resampled efficient frontier with shrinkage
    hold on; plot(vola_RE_shrink, rets_RE_shrink, 'c-', 'LineWidth',1.5)

    title('Efficient Frontiers')
    xlabel('Standard deviation $\sigma$', 'Interpreter','Latex')
    ylabel('Expected excess return $\mu^e$', 'Interpreter','Latex')
    grid on
    box on

    % assets
    hold on; plot(sqrt(diag(S_MLE)), mu, 'k*', 'MarkerSize',3)
    hold on; text(sqrt(diag(S_MLE))*1.01, mu, tickers)
    legend({'Mean-Variance MLE', 'Mean-Variance Shrinkage', ...
        sprintf('Resampled N=%i', N), sprintf('Resampled Shrinkage N=%i', N), 'Assets'}, ...
        'Location','southwest')


    % save to PDF
    set(gcf, 'Units','inches')
    pos = get(gcf, 'Position');   
    set(gcf, 'PaperUnits','inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
    print(gcf, '-dpdf', sprintf('../report/3_EFs_N=%i.pdf', N));


    
    % plot frontier portfolios composition
    figure
    set(gcf, 'Position',  [600, 400, 1000, 800])
    t = tiledlayout(2, 2, 'Padding','compact', 'TileSpacing','Compact');

    title(t, 'Frontier portfolios composition', 'FontSize',14, 'FontWeight','bold')
    xlabel(t, 'Rank $Z$', 'Interpreter','Latex')
    ylabel(t, 'Cumulative weight $\%$', 'Interpreter','Latex')

    nexttile
    area(1:Z, ef_MV_MLE.Weights*100)
    title('Mean-Variance MLE')
    axis([1, Z, 0, 100]);
    box on

    nexttile
    area(1:Z, ef_MV_shrink.Weights*100)
    legend(tickers, 'Location', 'NorthEastOutside')
    title('Mean-Variance Shrinkage')
    axis([1, Z, 0, 100]);
    box on

    nexttile
    area(1:Z, avg_wgts_RE*100)
    title(sprintf('Resampled N=%i', N))
    axis([1, Z, 0, 100]);
    box on

    nexttile
    area(1:Z, avg_wgts_RE_shrink*100)
    title(sprintf('Resampled Shrinkage N=%i', N))
    axis([1, Z, 0, 100]);
    box on


    % save to PDF
    set(gcf, 'Units','inches')
    pos = get(gcf, 'Position');   
    set(gcf, 'PaperUnits','inches');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
    print(gcf, '-dpdf', sprintf('../report/3_compositions_N=%i.pdf', N));

end





