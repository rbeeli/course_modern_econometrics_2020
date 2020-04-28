% ---------------------------------------------------
% Generates ACF and PACF plots of all tested models.
% ---------------------------------------------------

clear; clc; close all force; rng default;
addpath('lib')


% AR(4) model coefficients
a_n = 4; a = [repmat([0.4 -0.3 0.2], a_n, 1) linspace(-0.1, -0.6, a_n)'];

% MA(3) model coefficients
b_n = 3; b = [repmat(-0.1368, b_n, 1) linspace(-0.8673, -0.1956, b_n)' repmat(0.0046, b_n, 1)];





T = 2000;

% analyze AR(4) models
for i=1:a_n
    % plots
    if i == 1
        figure, set(gcf, 'Position',[300, 150, 800, 1000])
        t = tiledlayout(a_n, 2, 'Padding','compact');
        title(t, sprintf('AR(4)  a=[%.1f, %.1f, %.1f, a_4]', a(i, 1), a(i, 2), a(i, 3)))
        xlabel(t, 'lag')
    end
    
    % simulate AR(4)
    ts = armasim(T, 1, a(i, :));
    
    nexttile; autocorr(ts);
    title(sprintf('ACF  a_4 = %.1f', a(i, 4)));
    xlabel([]); ylabel([]);
    
    nexttile; parcorr(ts);
    title(sprintf('PACF  a_4 = %.1f', a(i, 4)));
    xlabel([]); ylabel([]);
end

Plot2Pdf('plots/models_AR(4).pdf')





% analyze MA(3) models
for i=1:b_n
    rr = MAroots(b(i, :));
    disp([abs(rr) NaN min(abs(rr))])
    
    % plots
    if i == 1
        figure, set(gcf, 'Position',[1200, 150, 800, 1000])
        t = tiledlayout(b_n, 2, 'Padding','compact');
        title(t, sprintf('MA(3)  b=[%.4f, b_2, %.4f]', b(i, 1), b(i, 3)))
        xlabel(t, 'lag')
    end
    
    % simulate MA(3)
    ts = armasim(T, 1, [], b(i, :));
    
    nexttile; autocorr(ts);
    title(sprintf('ACF  b_2 = %.3f  (rr_{min} = %.5f)', b(i, 2), min(abs(rr))))
    xlabel([]); ylabel([]);
    
    nexttile; parcorr(ts);
    title(sprintf('PACF  b_2 = %.3f  (rr_{min} = %.5f)', b(i, 2), min(abs(rr))))
    text(0.5, 0.9, sprintf('rr = [%.3f, %.3f, %.3f]', rr(1), rr(2), rr(3)))
    xlabel([]); ylabel([]);
end

Plot2Pdf('plots/models_MA(3).pdf')




% MA(3) to AR(20) 
figure
hold on
plot(cell2mat(arma2ar({}, b(1, :), 20)), 'r')
plot(cell2mat(arma2ar({}, b(2, :), 20)), 'b')
plot(cell2mat(arma2ar({}, b(3, :), 20)), 'g')
title('Coefficients of MA(3) converted to AR(20) model')
legend({strcat('b = [', join(string(b(1, :)), ', '), ']'), ...
    strcat('b = [', join(string(b(2, :)), ', '), ']'), ...
    strcat('b = [', join(string(b(3, :)), ', '), ']')}, ...
    'Location', 'southeast')

Plot2Pdf('plots/models_MA(3)_as_AR(20).pdf')




function rr = MAroots(b)
    % footnote on page 300 of Linear models and time-series
    % analysis book by Paolella, Marc S. (2019)
    rr = roots([b(end:-1:1) 1])';
end



