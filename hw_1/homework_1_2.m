clc; clear; close all force; rng default;


% parameters
n_loops = 5000000;                 % number of simulations to run / loops
p = [0.4; 0.6];                    % probability of having a daugther
r = [2; 3];                        % number of required daughters

x_min = sum(r);
x_max = 25;                        % max value for x-axis
x_vals = (x_min:x_max)';           % x axis values - total number of kids




% calculations

% analytic solution
analytic_vals_pmf = ana_pmf_vec(x_vals, p, r);
analytic_vals_cdf = ana_cdf_vec(x_vals, p, r);

% simulation
sim_hist = run_sim(n_loops, p, r, x_max); % returns histogram
sim_vals_pmf = sim_pmf_vec(x_vals, sim_hist);
sim_vals_cdf = sim_cdf_vec(x_vals, sim_hist);

% 1st order saddlepoint approximation
spa1_vals_pmf = spa1_pmf_vec(x_vals, p, r);
spa1_vals_cdf = spa1_cdf_vec(x_vals, p, r);

% 2nd order saddlepoint approximation
spa2_vals_pmf = spa2_pmf_vec(x_vals, p, r);





% plots
h = figure;
t = tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
title(t, sprintf('PMF and CDF estimation with r=%s and p=%s', ...
    strjoin(cellstr(num2str(r)), ','), ...
    strjoin(cellstr(num2str(p)), ',')));

% PMF
nexttile
plot(x_vals, sim_vals_pmf, 'Color','blue', 'LineWidth',1)
hold on
plot(x_vals, analytic_vals_pmf, 'Color','cyan', 'LineWidth',2, 'LineStyle', '-.')
hold on
plot(x_vals, spa1_vals_pmf, 'Color','magenta', 'LineWidth',2, 'LineStyle', '--')
hold on
plot(x_vals, spa2_vals_pmf, 'Color','green', 'LineWidth',2, 'LineStyle', ':')
title('Probability mass function (PMF)')
xlabel('$z$: Total number of kids', 'Interpreter', 'Latex')
ylabel('Probability: $f_Z(z) = P(Z = z)$', 'Interpreter', 'Latex')
legend({sprintf('Simulation (%i samples)', n_loops), 'Analytic', 'SPA 1st order', 'SPA 2nd order'}, 'Location', 'northeast')

% PMF relative percentage error
nexttile
plot(x_vals, pct_error(analytic_vals_pmf, sim_vals_pmf), 'Color','blue', 'LineWidth',1)
hold on
plot(x_vals, pct_error(analytic_vals_pmf, spa1_vals_pmf), 'Color','magenta', 'LineWidth',2)
hold on
plot(x_vals, pct_error(analytic_vals_pmf, spa2_vals_pmf), 'Color','green', 'LineWidth',2)
title('Probability mass function (PMF) - Relative percentage error')
ylim([-5 10])
xlabel('$z$: Total number of kids', 'Interpreter', 'Latex')
ylabel('Relative percentage error [%]')
legend({sprintf('Simulation (%i samples)', n_loops), 'SPA 1st order', 'SPA 2nd order'}, 'Location', 'northwest')

% CDF
nexttile
stairs(x_vals, sim_vals_cdf, 'Color','blue', 'LineWidth',1)
hold on
stairs(x_vals, analytic_vals_cdf, 'Color','cyan', 'LineWidth',2, 'LineStyle', '-.')
hold on
stairs(x_vals, spa1_vals_cdf, 'Color','magenta', 'LineWidth',1, 'LineStyle', '--')
title('Cumulative distribution function (CDF)')
ylim([-0.01 1.01])
xlabel('$z$: Total number of kids', 'Interpreter', 'Latex')
ylabel('$F_Z(z) = P(Z \leq z)$', 'Interpreter', 'Latex')
legend({sprintf('Simulation (%i samples)', n_loops), 'Analytic', 'SPA 1st order'}, 'Location', 'southeast')

% CDF relative percentage error
nexttile
plot(x_vals, pct_error(analytic_vals_cdf, sim_vals_cdf), 'Color','blue', 'LineWidth',1)
hold on
plot(x_vals, pct_error(analytic_vals_cdf, spa1_vals_cdf), 'Color','magenta', 'LineWidth',2)
title('Cumulative distribution function (CDF) - Relative percentage error')
ylim([-5 20])
xlabel('$z$: Total number of kids', 'Interpreter', 'Latex')
ylabel('Relative percentage error [%]')
legend({sprintf('Simulation (%i samples)', n_loops), 'SPA 1st order'}, 'Location', 'northwest')

% position on screen
set(gcf, 'Position',  [600, 200, 1400, 1000])

% save to pdf
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,'report/2_plots.pdf','-dpdf','-r0')


% print calculation results
headers  = {'Analytic'            'Simulation'           '1st order SPA'           '2nd order SPA'          };
res_ana  = [ana_pmf(8, p, r)      ana_pmf(12, p, r)      1 - ana_cdf(8, p, r)      1 - ana_cdf(12, p, r)    ];
res_sim  = [sim_pmf(sim_hist, 8)  sim_pmf(sim_hist, 12)  1 - sim_cdf(sim_hist, 8)  1 - sim_cdf(sim_hist, 12)];
res_spa1 = [spa1_pmf(8, p, r)     spa1_pmf(12, p, r)     1 - spa1_cdf(8, p, r)     1 - spa1_cdf(12, p, r)   ];
res_spa2 = [spa2_pmf(8, p, r)     spa2_pmf(12, p, r)     nan                       nan                      ];

err_sim  = pct_error(res_ana, res_sim);
err_spa1 = pct_error(res_ana, res_spa1);
err_spa2 = pct_error(res_ana, res_spa2);

% probabilities
res = [headers; [num2cell(res_ana') num2cell(res_sim') num2cell(res_spa1') num2cell(res_spa2')]];
res = [{''; 'P(Z=8)'; 'P(Z=12)'; 'P(Z>8)'; 'P(Z>12)'} res];

% relative percentage errors
res2 = [headers; [num2cell(zeros(size(res_ana))') num2cell(err_sim') num2cell(err_spa1') num2cell(err_spa2')]];
res2 = [{''; 'P(Z=8)'; 'P(Z=12)'; 'P(Z>8)'; 'P(Z>12)'} res2];

% pretty print tables
disp('Probabilities')
disp(' ')
disp(regexprep(evalc("disp(res)"), {'''','[',']'}, {'','',''}))

disp('Relative percentage errors')
disp(' ')
disp(regexprep(evalc("disp(res2)"), {'''','[',']'}, {'','',''}))






function errs = pct_error(true, approx)
    % Calculates the relative percentage error between
    % the true probability values and the approximated ones
    errs = 100*(approx - true) ./ bsxfun(@min, true, 1-true);
end

% --------------------------------------------------------------------------------------
% Simulation
% --------------------------------------------------------------------------------------

function hist = run_sim(n_loops, p, r, hist_max)
    % Simulates the setting by sampling from a uniform
    % distribution and counts the number of required trials
    % until the required number of daughters is reached for all r
    n_families = length(p);
    hist = zeros(hist_max, 1);

    for i=1:n_loops
       i_n_girls = zeros(n_families, 1);
       i_n_trials = zeros(n_families, 1);

       % loop until both have desired number of girls
       while any(i_n_girls < r)
          has_girl = rand(n_families, 1) <= p;
          is_trying = i_n_girls < r;
          i_n_girls = i_n_girls + (has_girl .* is_trying);
          i_n_trials = i_n_trials + is_trying;
       end

       % sanity checks (debug only)
       %assert(all(i_n_girls == r));
       %assert(all(i_n_trials >= i_n_girls));

       % update histogram
       n_combined_trials = sum(i_n_trials);
       if n_combined_trials <= length(hist)
           hist(n_combined_trials) = hist(n_combined_trials) + 1;
       end
    end
    
    % calculate probabilities
    hist = hist ./ sum(hist);
end

function p = sim_pmf(sim_hist, x)
    % Calculates the PMF value P(X = x) based on the histogram from run_sim(.)
    p = sim_hist(x);
end

function p = sim_cdf(sim_hist, x)
    % Calculates the CDF value P(X <= x) based on the histogram from run_sim(.)
    p = sum(sim_hist(1:x));
end

% vectorized versions of PMF/CDF functions
function p = sim_pmf_vec(x_vec, sim_hist)
    p = arrayfun(@(x) sim_pmf(sim_hist, x), x_vec);
end

function p = sim_cdf_vec(x_vec, sim_hist)
    p = arrayfun(@(x) sim_cdf(sim_hist, x), x_vec);
end



% --------------------------------------------------------------------------------------
% Analytic
% --------------------------------------------------------------------------------------

function p_res = ana_pmf(t, p, r)
    % Calculates the CDF value based on analytic probabilities (convolution).
    p_res = nan;
    
    % analytic solution only works for 2 r.v.s.
    if length(r) > 2
        return
    end
    
    % check argument is in support of r.v.
    if t < sum(r)
        return
    end
    
    p_res = 0;
    
    r1 = r(1);
    r2 = r(2);
    
    p1 = p(1);
    p2 = p(2);
    
    % sum over region for which t=x+y holds (convolution)
    for x=r1:(t-r2)
        y = t-x;
        p_x = nchoosek(x-1, r1-1) * p1^r1 * (1-p1)^(x-r1); % Negative Binomial
        p_y = nchoosek(y-1, r2-1) * p2^r2 * (1-p2)^(y-r2); % Negative Binomial
        p_res = p_res + p_x * p_y;
    end    
end

function p_res = ana_cdf(t, p, r)
    % Calculates the CDF value based on analytic probabilities (convolution).
    p_res = nan;
    
    % analytic solution only works for 2 r.v.s.
    if length(r) > 2
        return
    end
    
    p_res = 0;
    
    for i=sum(r):t
        p_res = p_res + ana_pmf(i, p, r);
    end
end

% vectorized versions of PMF/CDF functions
function p = ana_pmf_vec(x_vec, p, r)
    p = arrayfun(@(x) ana_pmf(x, p, r), x_vec);
end

function p = ana_cdf_vec(x_vec, p, r)
    p = arrayfun(@(x) ana_cdf(x, p, r), x_vec);
end


% --------------------------------------------------------------------------------------
% Saddlepoint approximation (SPA)
% --------------------------------------------------------------------------------------
    
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
    k = sum(-(((p-1)*exp(t) .* r .* ((p-1).^2*exp(2*t) - 4*(p-1)*exp(t) + 1)) ./ ((p-1)*exp(t) + 1).^4));
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


% Tried 2nd order S.P.A. for CDF according to Slides/Book, but the result is not stable / accurate,
% probably because of some singularities - further investiation required.

% function p_res = spa2_cdf(x, p, r)
%     % Second-order saddlepoint approximation of CDF
%
%     % Add +1 to x to account for strict inequality in SPA CDF formula,
%     % this obviously only works in the discrete case
%     x = x+1;
%     
%     p_res = nan;
%     
%     % check argument is in support of r.v.
%     if x <= sum(r)
%         return
%     end
%     
%     % 1st order SPA formula for CDF
%     [p_1st, sHat, w, u] = spa1_cdf(x, p, r);
%     
%     % 2nd order SPA formula for CDF
%     kappa3 = K3(sHat, p, r)/(K2(sHat, p, r)^(3/2));
%     kappa4 = K4(sHat, p, r)/(K2(sHat, p, r)^(4/2));
%     p_res = p_1st - normpdf(w)*( ...
%           u^-1*( 1/8*kappa4 - 5/24*kappa3^2 ) ...
%         - u^-3 ...
%         - kappa3/(2*u^2) ...
%         + 1/w^3 ...
%     );
% end


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







