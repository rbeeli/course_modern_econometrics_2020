function [R, dates, T] = LoadDJIARets()
    % Loads the file DJIA_daily.csv and
    % computes the percentage log-returns based on `Adj Close`.
    tbl = readtable('DJIA_daily.csv', 'PreserveVariableNames',true);
    dates = tbl.Date(2:end);
    
    % percentage log-returns
    lP = log(tbl.('Adj Close'));
    R = 100*(lP(2:end) - lP(1:(end-1)));
    T = length(R);
end