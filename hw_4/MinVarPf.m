function minVarPf = MinVarPf(expRets, covMat, allowShorts)
    % computes the minimum variance portfolio
    nAssets = length(expRets);
    x0 = (1/nAssets)*ones(nAssets, 1);
    Aeq = ones(1, nAssets);
    beq = 1;
    
    % weight bounds
    ub = ones(nAssets, 1);
    lb = zeros(nAssets, 1);
    if allowShorts
        lb = -ones(nAssets, 1);
    end
    
    % optimizer
    options = optimset('Algorithm','interior-point-convex','Display','off');
    wgts = quadprog(covMat, x0, [], [], Aeq, beq, lb, ub, [], options);
    
    minVarPf.Return = wgts'*expRets;
    minVarPf.Weights = wgts;
    minVarPf.Risk = sqrt(wgts'*covMat*wgts);
end
