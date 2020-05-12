function [var, var_violations, var_violations_pct, params, mixNsigma2] = MixN32GARCH11(R, w, var_lvls)
    % Estimate a Gaussian mixture model MixN(3,2) combined with GARCH(1,1)
    % variances using a moving window of size w.
    % Also compute the Value-at-Risk for the given significance
    % levels var_lvls.
    %
    % Parameters:
    %    R              Returns vector
    %    w              Moving window size
    %    var_lvls       Value-at-Risk significance levels
    
    disp('Parallel estimation of MixN(3,2)-GARCH(1,1)...')

    T = length(R);
    
    % estimated parameters, see GetMixNParams(x) for usage
    params = NaN(T, 12);
    
    % predicted variances of Gaussian Mixture model
    mixNsigma2 = NaN(T, 3);
    
    % mean return of each window
    wnd_means = NaN(T, 1);
    
    % number of failed fmincon optimizations
    nfailed = 0;
    
    parfor i=w+1:T
        % window data
        wnd = R((i - w):(i-1));
        
        % MixN-GARCH model was constructed to
        % have zero-mean, hence we subtract the
        % mean return of the estimation window
        % to get a mean-zero time series.
        % Store it for later VaR calculation.
        wnd_means(i) = mean(wnd);
        wnd = wnd - wnd_means(i);

        % estimate MixN(3,2)-GARCH(1,1)
        [m_params, ll, exitflag] = EstimateMixN32GARCH11(wnd);

        if exitflag == 2 || exitflag == 1
            % optimizer finished, store new parameters
            params(i, :) = m_params;
        else
            % optimizer failed, ignore parameters
            nfailed = nfailed + 1;
            fprintf('Optimizer exited with code %i at i=%i \n', exitflag, i)
        end

        fprintf('%i\n', i)
    end

    % forward-fill failed parameter estimations
    % use previous period parameters if optimizer fails (no look-ahead bias)
    fprintf('Optimizer failed %i out of %i times (%.2f%%) \n', nfailed, T-w, 100*nfailed/(T-w));
    params = fillmissing(params, 'previous');

    % predict volatilities and VaR
    disp('VaR prediction...')
    var = zeros(T, length(var_lvls));  % VaR

    parfor i=w+1:T
        % window data
        wnd = R((i - w):(i-1));

        % unpack weights from fmincon params vector
        [wgts, means, alphas, betas, gammas] = GetMixNParams(params(i,:));

        % predict volatility
        [~, sigma2pred] = Garch11sigma2(wnd, alphas, betas, gammas);
        mixNsigma2(i, :) = sigma2pred;

        % compute Value-at-Risk
        var(i,:) = MixNQuantile(var_lvls', wgts, means, sqrt(sigma2pred));
        
        % add back mean-return of estimation window for VaR
        var(i,:) = wnd_means(i) + var(i,:);
    end
    
    % VaR violations per VaR level
    var_violations = R*ones(1,3) < var;
    var_violations_pct = mean(var_violations(w+1:end, :));
    
    disp('Model estimation finished.')
end

function [params, loglik, exitflag] = EstimateMixN32GARCH11(R)
    % MixN(3,2)-GARCH(1,1) diagonal model.
    %
    % Parameters:
    %    R          Percentage log-returns
    %
    % Output:
    %    params     Vector with estimated parameters. Use
    %               GetMixNParams(params) for extraction.
    %    loglik     Log-likelihood of esimated parameters on R
    %    exitflag   Exit flag of fmincon optimizer
    
    x0 = [
        0.25, 0.70, 0.05, ...    % MixN weights
        0.1, 0.02, ...           % MixN means
        0.00001, 0.01, 1.8 ...   % GARCH alphas (constant)
        0.01, 0.16, ...          % GARCH betas (ARCH coef)
        0.9, 0.8, ...            % GARCH gammas (GARCH coef)
    ];
    opts = optimset('Display','notify', 'LargeScale','off', ...
                    'Maxiter',2000, 'TolFun',1e-6, 'TolX',1e-6);
    
    % lower bounds
    lb = [
        0, 0, 0, ...             % MixN weights
        -Inf, -Inf, ...          % MixN means
        0, 0, 0, ...             % GARCH alphas (constant)
        0, 0, ...                % GARCH betas (ARCH coef), 3rd=0
        0, 0, ...                % GARCH gammas (GARCH coef), 3rd=0
    ];

    % equality constraints
    Aeq = [1 1 1 0 0 0 0 0 0 0 0 0];
    beq = 1; % sum weights to 1
    
    % maximize log-likelihood using constrained optimizer
    fun = @(x) -MixNloglik(R, x);
    [params, loglik, exitflag] = fmincon(fun, x0, [], [], Aeq, beq, lb, [], [], opts);
end

function loglik = MixNloglik(R, x)
    T = length(R);
    
    % extract parameters
    [wgts, means, alphas, betas, gammas] = GetMixNParams(x);
    
    % compute GARCH variances
    sigma2 = Garch11sigma2(R, alphas, betas, gammas);
    
    % compute likelihoods
    likelihoods = zeros(T,3);
    for k=1:3
        likelihoods(:,k) = normpdf(R, ones(T,1)*means(k), sqrt(sigma2(:, k)));
    end
    likelihoods = likelihoods * wgts;
    
    % log-likelihood
    loglik = sum(log(likelihoods));
end

function p = MixNCDF(x, wgts, mu, sigma)
    % Cumulative distribution function
    % of Gaussian mixture model.
    % 
    % Parameters:
    %    x      Vector of values for F(x)
    %    wgts   Weights of Gaussian mixture dist.
    %    mu     Means of Gaussian mixture dist.
    %    sigma  Std. dev. of Gaussian mixture dist.
    
    p = wgts' * normcdf(x, mu, sigma);
end

function z = MixNQuantile(p, wgts, means, sigma)
    % Computes the quantile of a Gaussian mixture
    % distribution using a numerical optimizer.
    %
    % Parameters:
    %    p      Quantile to evaluate
    %    wgts   Weights of Gaussian mixture dist.
    %    mu     Means of Gaussian mixture dist.
    %    sigma  Std. dev. of Gaussian mixture dist.
    
    z = zeros(length(p),1);
    opts = optimset('Display','notify', 'TolX',1e-5);
    for i=1:length(p)
        fun = @(x) MixNCDF(x, wgts, means, sigma) - p(i);
        z(i) = fzero(fun, 0, opts);
    end
end

function [sigma2, sigma2pred] = Garch11sigma2(R, alphas, betas, gammas)
    % Computes the predicted sigma^2 for whole returns time series R.
    %
    % Parameters:
    %    alphas         Constant terms vector
    %    betas          ARCH terms vector
    %    gamms          GARCH terms vector
    %
    % Output:
    %    sigma2         GARCH sigma2 vector for 2:T
    %    sigma2pred     GARCH predicted sigma2 for T+1
    
    T = length(R);
    R2 = R.^2;
    sigma2 = ones(T, 1)*alphas';
    for t=2:T
        sigma2(t, :) = alphas + betas*R2(t-1) + gammas.*sigma2(t-1, :)';
    end
    
    sigma2pred = alphas + betas*R2(end) + gammas.*sigma2(end, :)';
end