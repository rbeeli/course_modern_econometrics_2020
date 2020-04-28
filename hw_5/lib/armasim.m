function y = armasim(T, sig2, pv, qv)

if nargin < 4, qv = []; end

p = length(pv); q = length(qv); pv = reshape(pv,1,p); qv = reshape(qv,1,q);
warmup = 1000; e = sqrt(sig2)*randn(T + warmup, 1); init = 0;
evec = zeros(q,1); yvec = zeros(p,1); y = zeros(T + warmup, 1);

for i=1:T+warmup
    if p > 0, y(i) = y(i) + pv*yvec; end
    if q > 0, y(i) = y(i) + qv*evec; end
    y(i) = y(i) + e(i);
    if p > 1, yvec(2:p) = yvec(1:p-1); end, yvec(1) = y(i);
    if q > 1, evec(2:q) = evec(1:q-1); end, evec(1) = e(i);
end

y = y(warmup + 1:end);

end


