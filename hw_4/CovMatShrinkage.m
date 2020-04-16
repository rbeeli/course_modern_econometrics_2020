function shrinkage = CovMatShrinkage(stocks_ex, F)
    % Calculates the optimal shrinkage coefficient.
    % Calculations are consistent with M.Wolf's
    % implementation covCor.m and his paper (tested):
    %
    %   Honey, I Shrunk the Sample Covariance Matrix
    %    - Journal of Portfolio Management, 2004.
    %
    
    [n, p] = size(stocks_ex);
    
    % demeaned returns
    r_demean = stocks_ex - mean(stocks_ex);

    % MLE of covariance matrix
    S_MLE = cov(stocks_ex, 1); % MLE 1/n
    S_MLE_diag = diag(S_MLE);

    gamma = sum(sum((F - S_MLE).^2));
    
    pi_ij = 1/n*(r_demean.^2)'*(r_demean.^2) - 2/n*(r_demean'*r_demean).*S_MLE + S_MLE.^2;
    pi = sum(sum(pi_ij));
    
    % 1/(p*(p-1)) instead of 2/(p*(p-1)) because I sum over
    %   the upper and lower triangle of symmetric matrix.
    % -p to remove diagonal elements (1s)
    sigma_bar = 1/(p*(p-1)) * (sum(sum(S_MLE./sqrt(S_MLE_diag(:, ones(p,1)).*S_MLE_diag(:, ones(p,1))'))) - p); 

    rho = 0;
    for i=1:p
        % diagonal elements
        rho = rho + pi_ij(i, i);
        
        % off-diagonal elements
        for j=1:p
            if j ~= i
                rho = rho + sigma_bar/2 * ( ...
                    sqrt(S_MLE(j,j)/S_MLE(i,i))*theta_kij(i, i, j, r_demean, S_MLE, n) + ...
                    sqrt(S_MLE(i,i)/S_MLE(j,j))*theta_kij(j, i, j, r_demean, S_MLE, n) ...
                );
            end 
        end
    end
    
    kappa = (pi - rho)/gamma;
    
    % optimal shrinkage coefficient
    shrinkage = min(1, max(0, kappa / n));
end


function out = theta_kij(k, i, j, r_demean, S, n)
    out = 0;    
    for t=1:n
       out = out + ( r_demean(t, k).^2 - S(k,k) )*( r_demean(t, i)*r_demean(t, j) - S(i, j) ); 
    end
    out = 1/n * out;
end
