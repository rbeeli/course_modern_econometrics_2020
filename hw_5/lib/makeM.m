function M = makeM(X)
    % Computes the regression residual maker matrix M.
    [T, ~] = size(X);
    M = eye(T) - X*inv(X'*X)*X';
end
