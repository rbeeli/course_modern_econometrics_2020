function [var, var_violations, var_violations_pct, params] = GARCH11(R, w, var_lvls, use_matlab_garch)
    % Estimate a Gaussian GARCH(1,1) model forvariances using a
    % moving window of size w. Also compute the Value-at-Risk for
    % the given significance levels var_lvls.
    %
    % Parameters:
    %    R                   Returns vector
    %    w                   Moving window size
    %    var_lvls            Value-at-Risk significance levels
    %    use_matlab_garch    Flag for using Econometrics Toolbox GARCH(1,1)
    %                        estimator or the one based on Paolella's book (babygarch).
    
    T = length(R);
    
    var = zeros(T, length(var_lvls));     % VaR
    params = NaN(T,4);
    nfailed = 0;

    parfor i=w+1:T
        % window estimation
        wnd = R((i - w):(i-1));

        if use_matlab_garch
            % estimate model using Matlab's Econometrics toolbox
            mdl = garch('Offset',NaN, 'GARCHLags',1, 'ARCHLags',1);
            estMdl = estimate(mdl, wnd, 'Display','off');
            mdlparams = [estMdl.Offset estMdl.Constant estMdl.ARCH{1} estMdl.GARCH{1}];
            params(i,:) = mdlparams;

            % predict volatility 1-step ahead
            %sigma2pred = forecast(estMdl, 1, 'Y0',wnd); % gives same result as below
            [~, sigma2pred] = Garch11sigma2(wnd, mdlparams);
        else
            % param = [returns mean, constant, ARCH coef, GARCH coef]
            [mdlparams, ~, ~, ~, exitflag] = babygarch(wnd);
            if exitflag == 0
                nfailed = nfailed + 1;
            end
            params(i,:) = mdlparams;

            % predict volatility 1-step ahead
            [~, sigma2pred] = Garch11sigma2(wnd, mdlparams);
        end

        % compute VaR breaches
        var(i,:) = norminv(var_lvls, mean(wnd), sqrt(sigma2pred));

        if mod(i, 50) == 0
            fprintf('%i\n', i)
        end
    end

    fprintf('Optimizer failed %i out of %i times (%.2f%%) \n', nfailed, T-w, 100*nfailed/(T-w));

    % VaR breaches per alpha level
    var_violations = R*ones(1,3) < var;
    var_violations_pct = mean(var_violations(w+1:end,:));
end


function [sigma2, sigma2pred] = Garch11sigma2(R, params)
    % Computes the predicted sigma^2 for whole returns time series R.
    %
    % Parameters:
    %    params         Vector of [mean return, constant, ARCH coef, GARCH coef]
    %
    % Output:
    %    sigma2         GARCH sigma2 vector for 2:T
    %    sigma2pred     GARCH predicted sigma2 for T+1
    
    T = length(R);
    e = R - params(1);
    e2 = e.^2;
    
    % compute GARCH variances sigma2
    sigma2 = ones(T, 1);
    sigma2(1) = params(2) + params(3)*mean(e2) + params(4)*mean(e2);
    for t=2:T
        sigma2(t) = params(2) + params(3)*e2(t-1) + params(4)*sigma2(t-1);
    end
    
    sigma2pred = params(2) + params(3)*e2(end) + params(4)*sigma2(end);
end