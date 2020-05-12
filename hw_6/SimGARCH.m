function [R, dates, T] = SimGARCH(T)
    % Simulates a GARCH returns time series
    % with heteroscedastic volatilities.
    % Based on book "Linear Models and Time-Series Analysis"
    % of Marc S. Paolella .
    
    % arbitrary dates from 2000-1-1 on
    dates = (datetime(2000,1,1):(datetime(2000,1,1) + days(T-1)))';
    
    a = randn(T+1,1)/12;
    e = 0.111 + randn(T+1,1)/1;

    P = zeros(T+1,1);
    P(1) = 100;
    for t=2:T
        P(t) = (1 + a(t))*P(t-1) + e(t);
    end

    if any(P <= 0.001)
        P = P + abs(min(P)) + 0.001;
    end
    
    % percentage log-returns
    lP = log(P);
    R = 100*(lP(2:end) - lP(1:(end-1)));
    
    % limit to max 50% movement
    R = max(-50*ones(T,1), min(R, 50*ones(T,1)));
end
