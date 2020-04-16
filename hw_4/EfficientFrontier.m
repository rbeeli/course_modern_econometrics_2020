function ef = EfficientFrontier(expRets, covMat, nPfs, allowShorts)
    % Calculates Mean-Variance efficient portfolios using the
    % supplied expected returns and covariance matrix. Also,
    % the number of portfolios on the frontier `nPfs` and whether
    % short-selling `allowShorts` is allowed, can be configured.
    
    % check that covariance matrix is positive definite,
    % otherwise quadprog won't be able to solve the problem
    [~, p] = chol(covMat);
    assert(p == 0);
    
    nAssets = length(expRets);

    % find minimum variance portfolio
    ef.MinVarPf = MinVarPf(expRets, covMat, allowShorts);

    % define expected returns range for efficient frontier
    maxRet = max(expRets);
    pfRets = linspace(ef.MinVarPf.Return, maxRet, nPfs)';
    pfVars = NaN(nPfs, 1);
    pfWgts = NaN(nPfs, nAssets);
    
    % weight bounds
    ub = ones(nAssets, 1);
    lb = zeros(nAssets, 1);
    if allowShorts
        lb = -ones(nAssets, 1);
    end
    
    % optimize for each rank/expected return
    options = optimset('Algorithm','interior-point-convex','Display','off');
    x0 = (1/nAssets)*ones(nAssets, 1);
    
    for i=1:nPfs
        expRet = pfRets(i);
        Aeq = [ones(1, nAssets); expRets'];
        beq = [1; expRet];
        wgts = quadprog(covMat, x0, [], [], Aeq, beq, lb, ub, [], options);
        pfWgts(i, :) = wgts;
        pfVars(i) = wgts'*covMat*wgts;
    end

    ef.Return = pfRets;
    ef.Risk = sqrt(pfVars);
    ef.Weights = pfWgts;
end
