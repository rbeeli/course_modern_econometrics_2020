%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Course:       Modern econometric and statistical learning
%               methods for quantitative asset management
%
% Instructor:   Prof. Dr. Marc Paolella, Urban Ulrych
%               University of Zurich
%
% Author:       Rino Beeli
%
% Date:         May 12th, 2020
% 
% Topic:        Homework 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all force; rng default;



% parameters
[R, dates, N] = LoadDJIARets();
Ts = [50, 100, 150, 200, 300, 500, 800, 1000, 1500, 2000, 3000, 4000];
samples = 200;



% estimation
coefs = zeros(length(Ts), samples);
parfor i=1:length(Ts)
    T = Ts(i);
    for j=1:samples
        idx = randi([1 N-T]); % random start point
        param = babygarch(R(idx:(idx+T-1)));
        coefs(i, j) = sum(param(3:4));
    end
end



% plots
figure
set(gcf, 'Position',[700, 200, 750, 1000])
tiledlayout(2,1,'Padding','compact')

nexttile, plot(Ts, mean(coefs, 2), 'k-')
title('Sum of GARCH coefficients')
xlabel('Sample size T')
ylabel('Sum of coefficients')
yline(1, 'r--');
ylim([0.5 1.1])
legend({'Sum GARCH coefs', 'IGARCH border'}, 'Location','southeast');

nexttile, plot(log(Ts), mean(coefs, 2), 'k-')
title('Sum of GARCH coefficients (log T)')
xlabel('log(Sample size T)')
ylabel('Sum of coefficients')
yline(1, 'r--');
ylim([0.5 1.1])
legend({'Sum GARCH coefs', 'IGARCH border'}, 'Location','southeast');

%Plot2Pdf('output/GARCH11_border.pdf')