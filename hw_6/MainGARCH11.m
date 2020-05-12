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
% Estimates a Gaussian GARCH(1,1) model for predicting variances
% using a moving window of size w in order to
% compute Value-at-Risk at levels 1%, 2.5% and 5%.
% Uses a simulated time series with heteroskedasticity
% see function SimGARCH.
% ----------------------


% % DJIA index returns
% [R, dates] = LoadDJIARets();
% name = 'DJIA';

% simulate percentage log-returns
[R, dates] = SimGARCH(10000);
name = 'Simulated TS';

w = 500;                              % window size
var_lvls = [0.01, 0.025, 0.05];       % VaR significance levels
use_matlab_garch = false;             % Matlab GARCH estimation or babygarch(y)



% returns parameters of GARCH(1,1) model
% [mu, c_0, c_1, d_1]
% starts from row w+1 on to have values
params = BlackboxGARCH11(R, dates, w, var_lvls, name, use_matlab_garch);