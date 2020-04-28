function [num_coefs, sim_coefs, sim_coefs_sig, sim_LR, sim_LRboot] = AR_SPA(X, Y, B, alpha_c, alpha_coefs)
    % Simulates the Saddlepoint method for estimating the number
    % of lags, p*, in an AR(p*) process.
    %
    % Parameters
    %  X                Design matrix, e.g. [1, t]
    %  Y                Response matrix, where each row corresponds
    %                   to a simulated time series.
    %  B                Number of bootstrap samples used for signed LRs
    %  alpha_c          Significance levels c used in ButPao(Y,X,c)
    %  alpha_coefs      Significance level for signed LR test.
    %                   Can be a vector of length p_max or a scalar.
    
    rep = size(Y, 1);
    p_max = length(alpha_c);
    
    % simulation results
    num_coefs = NaN(rep, 1);
    sim_coefs = NaN(rep, p_max);
    sim_coefs_sig = NaN(rep, p_max);
    sim_LR = NaN(rep, p_max);
    sim_LRboot = cell(rep, 1);

    % parallel processing for faster calculation
    parfor r=1:rep
        fprintf('SPA rep=%i of %i \n', r, rep)
        
        % fit AR(p*) using SPA method
        [p_star, coefs, coefs_sig, LR, LRboot] = FitAR_SPA(X, Y(r, :)', B, alpha_c, alpha_coefs);

        if p_star == 0
            disp('Warning: No AR(p) lag found'); 
        else
            % store LR, LRboot, p*, AR(p*) coefs., coefs. test
            sim_LRboot{r} = LRboot;
            num_coefs(r) = p_star;

            tmp = NaN(1, p_max); tmp(1:p_star) = LR;
            sim_LR(r, :) = tmp;

            tmp = NaN(1, p_max); tmp(1:p_star) = coefs;
            sim_coefs(r, :) = tmp;

            tmp = NaN(1, p_max); tmp(1:p_star) = coefs_sig;
            sim_coefs_sig(r, :) = tmp;
        end
    end
    
    % drop data where SPA didn't find any AR(p) lag (almost never happens)
    sim_coefs = sim_coefs(~isnan(num_coefs), :);
    sim_coefs_sig = sim_coefs_sig(~isnan(num_coefs), :);
    sim_LR = sim_LR(~isnan(num_coefs), :);
    sim_LRboot = sim_LRboot(~isnan(num_coefs), :);
    num_coefs = num_coefs(~isnan(num_coefs));
end



function coefs_sig = SignedLRtest(LRboot, LRsim, alpha)
    % Tests wheter the bootstrapped signed likelihood ratios
    % reject the null hypothesis (zero coefficient) or not.
    %
    % Parameters:
    %  LRboot       Bootstrapped signed likelihood ratios
    %  LRsim        Likelihood ratios using real data
    %  alpha        Significance level (scalar or vector)
    p = length(LRsim);
    
    % convert alpha to vector if scalar
    if length(alpha) == 1
        alpha = ones(1, p) * alpha;
    end
    
    coefs_sig = NaN(1, p);
    for i=1:p
        coefs_sig(i) = quantile(LRboot(:, i), 1-alpha(i)/2) < LRsim(i) || quantile(LRboot(:, i), alpha(i)/2) > LRsim(i);
    end
end



function [p_star, coefs, coefs_sig, LR, LRboot] = FitAR_SPA(X, Y, B, alpha_c, alpha_coefs)
    % Uses the SPA method from chapter 9.5 in order
    % to determine AR(p*) and additionally, it Bootstraps the
    % signed LR statistics which can be used to determine
    % which coefficients are deemed zero.
    %
    % Parameters:
    %  X            Design matrix (tested with constant + time trend)
    %  Y            Response variable, residuals considered to follow AR(p*)
    %  B            Number of bootstrap samples used for signed LRs
    %  p_max        Upper limit for p*. p-values are linearly interpolated
    %               using linspace(0.175, 0.025, p_max) for c in ButPao(Y, X, c)
    %
    % Returns:
    %  p_star       Optimal p* of AR(p*)
    %  coefs        Estimated AR(p*) coefficients
    %  LRsim        Signed LR of all coefficients using the provided data
    %  LRboot       All bootstrapped signed LRs of all coefficients
    %
    
    T = length(Y);
    p_max = length(alpha_c);
    
    % ensure alpha_coefs either scalar or same length as alpha_c
    assert(p_max == length(alpha_coefs) || length(alpha_coefs) == 1);
    
    % OLS residuals (w/ autocorrelation)
    M = makeM(X);
    e = M*Y;
    
    % optimal p* based on SPA
    [~, p_star] = ButPao(e, [], alpha_c);
    p_star = min(p_star, p_max); % set upper bound to p_max (shouldn't happen though)
    if p_star == 0
        % SPA method did not detect an AR(p) process
        LR = []; LRboot = []; coefs = []; coefs_sig = [];
        return;
    end
    
    % fit AR(p*)
    [param_unr, ~, resid_unr, ~, ll_unr] = armareg(e, [], p_star, 0, 1);
    coefs = param_unr(1:p_star);
    
    % signed LR statistics for every i-th coef. set to 0
    LR = zeros(p_star, 1);
    for i=1:p_star
        ll_res = armareg_restricted(e, [], p_star, i, 0, 1);
        LR(i) = sign(coefs(i)) * sqrt(-2*(ll_res - ll_unr));
    end
    
    % bootstrap signed LR statistics for every i-th coef. set to 0
    LRboot = zeros(B, p_star);
    for b=1:B
        % non-parametric bootstrap
        resid_b = resid_unr(randi([1 T], T, 1));
        e_b = e - resid_unr + resid_b;
        
        % unrestricted MLE of AR(p*)
        [B_param_unr, ~, ~, ~, B_ll_unr] = armareg(e_b, [], p_star, 0, 1);
        B_coefs = B_param_unr(1:p_star);

        % signed likelihood ratio statistic for every i-th coef. set to 0
        for i=1:p_star
            B_ll_res = armareg_restricted(e_b, [], p_star, i, 0, 1);
            LRboot(b, i) = sign(B_coefs(i))*sqrt(-2*(B_ll_res - B_ll_unr));
        end
    end
    
    % test coefficients=0 based on signed LR statistics
    coefs_sig = SignedLRtest(LRboot, LR, alpha_coefs);
end



