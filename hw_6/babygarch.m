% --------------------------------------------------------------------------
% Codes for estimating a GARCH(1,1) model, based on book
% "Linear Models and Time-Series Analysis"
% of Marc S. Paolella with modified fminunc solver options.
% --------------------------------------------------------------------------

function [param, stderr, loglik, zvec, exitflag] = babygarch(y)
    % Normal-GARCH(1,1) with power=2
    % y is vector of log percentage returns
    initvec = [0 0.04 0.1 0.8];
    
    % mu c_0 c_1 d_1
    bound.lo    = [-4 0   0 0 ];
    bound.hi    = [ 4 0.5 1 1 ];
    bound.which = [ 1 1   1 1 ];
    opt = optimset('Display','notify', 'LargeScale','off', ...
                   'MaxFunEvals',10000, 'Maxiter',5000, ...
                   'TolFun',1e-7, 'TolX',1e-7);
    init = einschrk(initvec, bound);
    [pout,~,exitflag,~,~,hess] = fminunc(@(param) like(param,y,bound), init, opt);
    [loglik,zvec] = like(pout,y,bound);
    V = pinv(hess)/length(y);
    [param,V] = einschrk(pout,bound,V);
    stderr = sqrt(diag(V));
end

function [loglik,zvec] = like(param,y,bound)
    param = einschrk(real(param), bound, 999);
    meanterm = param(1);
    c0 = param(2);
    c1 = param(3);
    d1 = param(4);
    e = y - meanterm;
    [zvec, sigvec] = ungarch(e,c0,c1,d1);
    K = sqrt(2*pi);
    ll = -0.5 * zvec.^2 - log(K) - log(sigvec);
    loglik = -mean(ll);
end

function [eout, sigvec] = ungarch(e,c0,c1,d1)
    sigvec = zeros(length(e),1);
    e2 = e.^2;
    denom = 1-c1-d1;

    if denom > 0.001
        sinit = c0/denom;
    else
        sinit = mean(e2);
    end

    einit = sinit;

    % do the recursion in sigvec delta because it is faster
    sigvec(1) = c0+c1*einit+d1*sinit;
    for t=2:length(e)
        sigvec(t) = c0 + c1 *e2(t-1) + d1*sigvec(t-1);
    end

    sigvec = sigvec.^(1/2);
    eout = e./sigvec;
end

function [pout, Vout] = einschrk(pin, bound, Vin)
    lo = bound.lo;
    hi = bound.hi;
    welche = bound.which;

    if nargin < 3
        trans = sqrt((hi - pin) ./ (pin - lo));
        pout = (1 - welche) .* pin + welche .* trans;
        Vout = [];
    else
        trans = (hi + lo.*pin.^2) ./ (1 + pin.^2);
        pout = (1-welche) .* pin + welche .* trans;

        % now adjust the standard errors
        trans = 2 * pin .* (lo - hi) ./ (1 + pin.^2).^2;
        d = (1-welche) + welche .* trans; % either unity or delta method.
        J = diag(d);
        Vout = J * Vin * J ;
    end
end