function PlotVaR(titleStr, R, dates, w, var, var_lvls, var_viol, var_viol_pct)
    % Plots the returns time series for each VaR level including the
    % estimated VaRs (blue) and possible violations (red dots).
    % The percentage violations per VaR level is added as text.
    
    T = length(R);
    
    figure
    set(gcf, 'Position',[700, 200, 900, 1100])
    t = tiledlayout(length(var_lvls), 1,'Padding','compact');
    title(t, titleStr);
    
    for i=1:length(var_lvls)
        nexttile
        plot(dates, R, 'k-')
        title(sprintf('%.1f%% VaR', 100*var_lvls(i)))
        xlim([dates(w+1) dates(T)])
        ylim([1.1*min(min(var)) 1.3*max(R)])
        
        % VaR line
        hold on, plot(dates, var(:,i), 'b-')
        
        % VaR violations
        hold on, plot(dates(var_viol(:, i)), min(R), 'r*', 'MarkerSize',3)
        
        % Pct. violations text
        hold on, text(dates(w+round(0.03*(T-w))), max(R), sprintf('VaR violations: %.2f%%', 100*var_viol_pct(i)))
    end
end
