function [num_coefs, coefs, coefs_sig] = AR_LASSO(X, Y, p_max)
    % Uses LASSO to fit an AR(p) model, where p is the number of
    % non-zero coefficients at the lambda which minimizes
    % the MSE using Cross-Validation.
    %
    % Parameters
    %  X                Design matrix, e.g. [1, t]
    %  Y                Response matrix, where each row corresponds
    %                   to a simulated time series.
    %  p_max            Maximum number of possible lags (regressors).

    rep = size(Y, 1);
    T = size(Y, 2);
    num_coefs = NaN(rep, 1);
    coefs = NaN(rep, p_max);
    coefs_sig = NaN(rep, p_max);
    
    parfor r=1:rep
        fprintf('LASSO rep=%i of %i \n', r, rep)

        % OLS residuals (w/ autocorrelation)
        M = makeM(X);
        e = M*Y(r, :)';

        % design matrix for fitting ( e_{t-1}, ... , e_{t-p_max} ) onto e_t
        X_e = LassoDesignMatrix(e, p_max);
        Y_e = e((p_max+1):end);

        if 1==1
            % k-fold Cross-Validation
            [l_coefs, FitInfo] = lasso(X_e, Y_e, 'CV',5, 'Lambda', exp(linspace(-5, 3, 500)));
            lambdas = FitInfo.Lambda;

            % get min MSE coefficients
            lambda_idx = FitInfo.IndexMinMSE;
            lasso_coefs = l_coefs(:, lambda_idx);
        else
            % 50:50 split approach
            % worse than CV, too little data
            T_e = T - p_max;
            frac_train = 0.5;
            n_train = round(frac_train * T_e);
            n_test = T_e - n_train;
            X_e_train = X_e(1:n_train, :);
            Y_e_train = Y_e(1:n_train, :);
            X_e_test = X_e((T_e-n_test+1):T_e, :);
            Y_e_test = Y_e((T_e-n_test+1):T_e, :);

            % fit on training data
            [l_coefs, FitInfo] = lasso(X_e_train, Y_e_train);
            lambdas = FitInfo.Lambda; nLambdas = length(lambdas);

            % calculate MSE for all lambdas on test data
            lambda_MSEs = NaN(nLambdas, 1);
            for i=1:nLambdas
                i_pred = X_e_test*l_coefs(:, i) + FitInfo.Intercept(i);
                lambda_MSEs(i) = mean((Y_e_test - i_pred).^2);
            end

            % get min MSE coefficients
            lambda_idx = find(lambda_MSEs == min(lambda_MSEs), 1, 'last');
            lasso_coefs = l_coefs(:, lambda_idx);
        end
        
        if lambda_idx == length(lambdas)
            disp('Warning: Minimum MSE lambda at boundary (no CV convergence)');
        else
            % remove zero coefficients
            lasso_coefs(abs(lasso_coefs) < 0.0001) = NaN;

            % store results
            coefs(r, :) = lasso_coefs;
            coefs_sig(r, :) = ~isnan(lasso_coefs);
            num_coefs(r) = sum(~isnan(lasso_coefs));
        end
    end
    
    % drop data where LASSO didn't find CV minimum (happens when sample is very small)
    coefs = coefs(~isnan(num_coefs), :);
    coefs_sig = coefs_sig(~isnan(num_coefs), :);
    num_coefs = num_coefs(~isnan(num_coefs));
end


function X = LassoDesignMatrix(resid, p_max)
    % Creates a design matrix based on residuals,
    % where the time dependency is kept.
    % Each row i of X has values resid_{i+p_max-1}, ..., resid_{i}.
    % This way, the serial dependence of p_max lags is preserved.
    T = length(resid);
    X = NaN(T-p_max, p_max);
    for i=1:T-p_max
        % t-1 residuals in first column, t-2 in second, ...
        X(i, :) = flip(resid(i:(i + p_max - 1))'); % reverse order
    end
end
