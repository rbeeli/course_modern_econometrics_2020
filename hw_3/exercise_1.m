%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course:       Modern econometric and statistical learning
%               methods forquantitative asset management
%
% Instructor:   Prof. Dr. Marc Paolella, Urban Ulrych
%               University of Zurich
%
% Author:       Rino Beeli
%
% Date:         April 7th, 2020
% 
% Topic:        Homework 3 - Exercise 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all force; rng default;


% simulate returns
p = 4;
ts = randn(50, p) + 0.001 * mean(randn(50,1));

Sigma = cov(ts);            % covariance matrix
mu = mean(ts)';             % expected returns
R_f = 0.01;                 % risk-free rate
mu_e = mu - R_f;            % expected excess returns

% tangency portfolio
w_t = (inv(Sigma)*mu_e)/(ones(4,1)'*inv(Sigma)*mu_e);
mu_t = w_t'*mu_e;
sigma_t = w_t'*Sigma*w_t;

% generate plot data
mu_p_target = (-0.2:0.005:0.2)';
mu_p_opt = NaN(length(mu_p_target), 1);
sigma_p = NaN(length(mu_p_target), 1);

for i=1:length(mu_p_target)
    mu_p = mu_p_target(i);
    w_p = ((mu_p - R_f)/(mu_e'*inv(Sigma)*mu_e)) * (inv(Sigma) * mu_e);
    mu_p_opt(i) = w_p'*mu + (1 - sum(w_p))*R_f;
    sigma_p(i) = w_p'*Sigma*w_p;
end


% simulate portfolios with random weights
N = 10000;
mu_p_rnd = NaN(N, 1);
sigma_p_rnd = NaN(N, 1);

for i=1:N
    w_p = rand(4, 1)*0.8-0.4;
    mu_p_rnd(i) = w_p'*mu + (1 - sum(w_p))*R_f;
    sigma_p_rnd(i) = w_p'*Sigma*w_p;
end


% plot
figure
set(gcf, 'Position',  [600, 400, 1400, 400])
t = tiledlayout(1, 3, 'Padding','compact');

nexttile
plot(mu_p_target, mu_p_opt, 'r', 'LineWidth',1.5)
title('Expected returns')
xlabel('Expected portfolio return $\bf{\mu}_p$', 'Interpreter','latex')
ylabel('Optimized expected portfolio return')
ylabel('Optimized expected portfolio return $E[\mu_p^*]$', 'Interpreter','latex')

nexttile
plot(mu_p_target, sigma_p, 'r', 'LineWidth',1.5)
hold on
plot(mu_p_rnd, sigma_p_rnd, 'bo', 'MarkerSize',1)
%hold on 
%plot([mu_t], [sigma_t], 'go', 'MarkerSize',3)
title('Variance')
xlabel('Expected portfolio return $\bf{\mu}_p$', 'Interpreter','latex')
ylabel('(Optimized) portfolio variance $\sigma_p^2$', 'Interpreter','latex')

nexttile
plot(mu_p_target, sqrt(sigma_p), 'r', 'LineWidth',1.5)
hold on
plot(mu_p_rnd, sqrt(sigma_p_rnd), 'bo', 'MarkerSize',1)
title('Standard deviation')
xlabel('Expected portfolio return $\bf{\mu}_p$', 'Interpreter','latex')
ylabel('Optimized portfolio standard deviation $\sigma_p$', 'Interpreter','latex')

% set(gcf, 'Units','inches')
% pos = get(gcf, 'Position');   
% set(gcf, 'PaperUnits','inches');
% set(gcf, 'PaperSize', [pos(3) pos(4)]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);
% print(gcf, '-dpdf', 'report/1_a_1.pdf');
