function params = BlackboxGARCH11(R, dates, w, var_lvls, name, use_matlab_garch)
    % Estimates GARCH(1,1)-based VaR and violations and creates plots.
    %
    % Parameters:
    %   R                   Vector of returns
    %   dates               Vector of dates
    %   w                   Moving window size
    %   var_lvls            Vector of Value-at-Risk significance levels
    %   name                Name of time series (used in plot)
    %   use_matlab_garch    Switch for whether to use babygarch(y)
    %                       or Matlab's garch/estimate functions.
    
    % estimation
    [var, var_viol, var_viol_pct, params] = ...
        GARCH11(R, w, var_lvls, use_matlab_garch);

    % plots
    PlotVaR(sprintf('%s - GARCH(1,1)', name), ...
        R, dates, w, var, var_lvls, var_viol, var_viol_pct)

    %Plot2Pdf(sprintf('output/%s_GARCH11_VaR.pdf', name))
end
