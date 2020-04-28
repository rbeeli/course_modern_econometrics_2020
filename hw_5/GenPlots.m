function GenPlots(method, filename, model_desc, num_coefs, coefs, coefs_sig, p_max, LRboot)
    % Generates nice plots and saves them to a PDF file. The parameters of
    % this function are the return values of AR_LASSO(.) and AR_SPA(.).
    % 
    % Parameters:
    %  method       Either 'LASSO' or 'SPA'
    %  path         Path to target folder for PDF files
    %  filename     Filename prefix for all plots
    %  model_desc   Model description string
    %  num_coefs    Vector with number of selected coefs per simulation
    %  coefs        Matrix of estimated coefs per simulation
    %  coefs_sig    Matrix indicating if estimated coef is significant
    %  p_max        Configured p_max value
    %  LRBoot       Matrix of estimated signed likelihood-ratios
    %               using bootstrapping (SPA only)
    
    avg_coefs = nansum(coefs)/length(num_coefs);

    figure, set(gcf, 'Position',[1100, 200, 400, 900])
    t = tiledlayout(4, 1, 'Padding','compact');
    title(t, model_desc, 'FontSize',11);


    nexttile
    hist_dat = tabulate(num_coefs);
    bar(hist_dat(:, 3), 0.5)
    title(sprintf('Histogram p^* (%s)', method))
    xlabel('Number of lags')
    ylabel('%')
    
    
    dat = zeros(p_max, 1);
    for i=1:p_max
        sig = nansum(coefs_sig, 2);
        dat(i) = mean(sig(num_coefs == i));
    end
    nexttile
    bar(dat, 0.5); hold on; plot(0:p_max, 0:p_max, 'k--')
    title(sprintf('Mean significant coefs (%s)', method))
    xlabel('Number of lags')
    ylim([0 p_max])


    nexttile
    bar(mean(coefs_sig == 1)*100, 0.5)
    title(sprintf('Fraction of times coefficient selected (%s)', method))
    xlabel('Coefficient #')
    ylabel('%')
    

    nexttile
    bar(avg_coefs, 0.5)
    title(sprintf('Average coefficient value (%s)', method))
    xlabel('Coefficient #')
    tmp = avg_coefs; tmp(avg_coefs <= 0) = nan;
    text(1:p_max, tmp , num2str(tmp', 2), 'vert','bottom','horiz','center', 'FontSize',8);
    tmp = avg_coefs; tmp(avg_coefs > 0) = nan;
    text(1:p_max, tmp , num2str(tmp', 2), 'vert','top','horiz','center', 'FontSize',8);

    Plot2Pdf(sprintf('plots/%s_coefs_stats.pdf', filename))
    
    
    
    
    figure, set(gcf, 'Position',[1200, 200, 400, 900])
    t = tiledlayout(ceil(p_max/2), 2, 'Padding','compact');
    title(t, model_desc, 'FontSize',11);
    
    for i=1:p_max
        nexttile, histogram(coefs(:, i), 30);
        title(sprintf('Distribution of a_%i', i))
    end
    Plot2Pdf(sprintf('plots/%s_coefs_dist.pdf', filename))
    
    
    
    if strcmp(method, 'SPA')
        % plots distribution of bootstrapped LRs
        LRboot_ = LRboot{1};
        p = size(LRboot_, 2);

        % QQ-plots
        figure, set(gcf, 'Position',[1200, 300, 600, 450])
        tiledlayout(ceil(p/2), 2, 'Padding','compact')
        for i=1:p
            nexttile, qqplot(LRboot_(:, i));
            title(sprintf('QQ plot - LR stat of a_{%i}', i))
            ylabel('Sample Quantiles')
        end
        Plot2Pdf(sprintf('plots/%s_LRT_QQplot.pdf', filename))
    end
end
