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
clear; clc; close all force; rng(8);

% ----------------------
% Estimates a Gaussian mixture model MixN(3,2) combined GARCH(1,1)
% variances using a moving window of size w in order to
% compute Value-at-Risk at levels 1%, 2.5% and 5%.
% Uses daily Dow Jones IA index returns from
% January, 2000 until May, 2020.
% ----------------------


[R, dates] = LoadDJIARets();    % DJIA index returns
name = 'DJIA';

% % simulate percentage log-returns
% [R, dates] = SimGARCH(10000);
% name = 'Simulated TS';

w = 500;                           % window size
var_lvls = [0.01, 0.025, 0.05];    % VaR significance levels



% returns parameters of MixN density
% [w_1, w_2, w_3, mu_1, mu_2, mu_2, sig2_1, sig2_2, sig2_3]
% starts from row w+1 on to have values
params = BlackboxMixN(R, dates, w, var_lvls, name);