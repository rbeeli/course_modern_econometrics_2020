clc; clear; close all force; rng default;


% parameters
p = [0.4; 0.6];             % probability of having a daugther
r = [2; 3];                 % number of required daughters
x = (sum(r):30)';           % x values


% 1st order saddlepoint approximation
spa1_vals_pmf = spa1_pmf_vec(x, p, r);
spa1_vals_cdf = spa1_cdf_vec(x, p, r);

% 2nd order saddlepoint approximation
spa2_vals_pmf = spa2_pmf_vec(x, p, r);



% cumulant generating function and derivatives 1 to 4
function k = K(t, p, r) 
    k = sum(r .* log(p) + r*t - r .* log(1 - (1-p)*exp(t)));
end

function k = K1(t, p, r)
    k = sum(r - ((p-1)*exp(t) .* r) ./ ((p-1)*exp(t) + 1));
end

function k = K2(t, p, r)
    k = sum(-(((p-1)*exp(t) .* r) ./ ((p-1)*exp(t) + 1).^2));
end

function k = K3(t, p, r)
    k = sum(((p-1)*exp(t) .* r .* ((p-1)*exp(t) - 1)) ./ ((p-1)*exp(t) + 1).^3);
end

function k = K4(t, p, r)
    k = sum(-((p-1)*exp(t) .* r .* ((p-1).^2*exp(2*t) - 4*(p-1)*exp(t) + 1)) ./ ((p-1)*exp(t) + 1).^4);
end



function sHat = spa_sHat(x, p, r)
    % Calculates sHat by solving s in the equation K'(s) - x = 0.
    % Around the convergence strip of the m.g.f., there
    % exists a unique root for the equation.
    
    % convergence strip
    cs = [-100, 0.999999*min(-log(1-p))];
    
    % solve K1[sHat] = x
    opts = optimset('Display','none', 'TolX',1e-16);
    func = @(s) K1(s, p, r) - x;
    sHat = fzero(func, cs, opts);
end


function [p_res, sHat] = spa1_pmf(x, p, r)
    % First-order saddlepoint approximation of PMF
    p_res = nan;
    sHat = nan;
    
    % check argument is in support of SPA for r.v.
    % we could use the SPA of the CDF to get a value for
    % the bounds, but we leave it like this.
    if x <= sum(r)
        return
    end
    
    % solve s in K'(s) = x
    sHat = spa_sHat(x, p, r);
    
    % 1st order SPA formula
    p_res = (1 / sqrt(2 * pi * K2(sHat, p, r))) * exp(K(sHat, p, r) - x * sHat);
end


function p_res = spa2_pmf(x, p, r)
    % Second-order saddlepoint approximation of PMF
    p_res = nan;
    
    % check argument is in support of r.v.
    if x <= sum(r)
        return
    end
    
    % 1st order SPA formula
    [p_1st, sHat] = spa1_pmf(x, p, r);
    
    % 2nd order SPA formula
    k2 = K2(sHat, p, r);
    k3 = K3(sHat, p, r);
    k4 = K4(sHat, p, r);
    p_res = p_1st * (1 + 1/8*( k4/(k2^2) ) - 5/24*( k3/(k2^(3/2)) )^2);
end


function [p_res, sHat, w, u] = spa1_cdf(x, p, r)
    % First-order saddlepoint approximation of CDF
    
    % Add +1 to x to account for strict inequality in SPA CDF formula,
    % this obviously only works in the discrete case!
    x1 = x+1;
    
    p_res = nan;
    sHat = nan;
    w = nan;
    u = nan;
    
    % check argument is in support of r.v.
    if x1 <= sum(r)
        return
    end
    
    % solve s in K'(s) = x
    sHat = spa_sHat(x1, p, r);
    
    % 1st order SPA formula for CDF
    w = sign(sHat) * sqrt(2*sHat*x1 - 2*K(sHat, p, r));
    
    % There exists a singularity when X=E[X], because w becomes 0.
    % Remove it via linear interpolation around X+/-0.1 when X is close to E[X]
    if abs(K1(0, p, r) - x1) < 1e-6
        p_res = 1/2*(spa1_cdf(x-0.1, p, r) + spa1_cdf(x+0.1, p, r));
    else
        u = (1 - exp(-sHat)) * sqrt(K2(sHat, p, r));
        p_res = normcdf(w) + normpdf(w)*(1/w - 1/u);
    end
end


% vectorized versions of PMF/CDF functions
function p = spa1_pmf_vec(x_vec, p, r)
    p = arrayfun(@(x) spa1_pmf(x, p, r), x_vec);
end

function p = spa1_cdf_vec(x_vec, p, r)
    p = arrayfun(@(x) spa1_cdf(x, p, r), x_vec);
end

function p = spa2_pmf_vec(x_vec, p, r)
    p = arrayfun(@(x) spa2_pmf(x, p, r), x_vec);
end