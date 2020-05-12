function params = BlackboxMixN(R, dates, w, var_lvls, name)
    % Estimates MixN(3,2)-GARCH(1,1)-based VaR and violations and creates plots.
    %
    % Parameters:
    %   R                   Vector of returns
    %   dates               Vector of dates
    %   w                   Moving window size
    %   var_lvls            Vector of Value-at-Risk significance levels
    %   name                Name of time series (used in plot)

    % estimation
    [var, var_viol, var_viol_pct, mparams, mixNsigma2] = ...
        MixN32GARCH11(R, w, var_lvls);

    % extract parameters
    T = length(R);
    params = NaN(T, 9);
    for i=w+1:T
        [wgts, means] = GetMixNParams(mparams(i,:));
        params(i,:) = [wgts' means' mixNsigma2(i,:)]; 
    end
    
    % plots
    PlotVaR(sprintf('%s - MixN(3,2)-GARCH(1,1)', name), ...
            R, dates, w, var, var_lvls, var_viol, var_viol_pct)

    %Plot2Pdf(sprintf('output/%s_MixN32-GARCH11_VaR.pdf', name))
end