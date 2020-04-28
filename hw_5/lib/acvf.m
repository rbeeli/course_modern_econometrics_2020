function Sigma=acvf(a,b,T);
% Sigma=acvf(a,b,T);
% a=(a_1,...,a_p) with sign convention: Y_t = a_1 Y_{t-1} + a_2  Y_{t-2} + ...
% b=(b_1,...,b_q), and T is desired size of the covariance matrix
% assumes T>max(p,q)

p=length(a); q=length(b); m=max(p,q)+1; %+1 because of the way leeuwARMA works
V=leeuwARMA(a,b,m); gamma=[V(:,1) ; zeros(T-m,1) ];
for k=m+1:T
  s=0; for i=1:p, s=s+a(i)*gamma(k-i); end, gamma(k)=s;
end
Sigma=toeplitz(gamma);

function V=leeuwARMA(a,b,T);
% a=(a_1,...,a_p) with sign convention: Y_t = a_1 Y_{t-1} + a_2  Y_{t-2} + ...
% b=(b_1,...,b_q), and T is desired size of the covariance matrix
% designed just for mixed ARMA processes. See leeuwAR and leeuwMA for the special cases.

p=length(a); q=length(b); a=reshape(a,1,p); b=reshape(b,1,q);

% Now zero pad
m=max(p,q); a=[a zeros(1,m-p)]; b=[b zeros(1,m-q)]; p=m; q=m;

a=-a; Tuse=T+p;
firrowP = [1 zeros(1,Tuse-1)]; fircolP = [1 a zeros(1,Tuse-p-1)];
P = toeplitz(fircolP,firrowP);
firrowQ1 = a(p:-1:1); fircolQ1 = [a(p) zeros(1,p-1)];
Q1 = toeplitz(fircolQ1,firrowQ1); Q = [Q1; zeros(Tuse-p,p)];

firrowM = [1 zeros(1,T-1)]; fircolM = [1 b zeros(1,T-q-1)];
M = toeplitz(fircolM,firrowM);
firrowN1 = b(q:-1:1); fircolN1 = [b(q) zeros(1,q-1)];
N1 = toeplitz(fircolN1,firrowN1); N = [N1; zeros(T-q,q)];

middle=P'*P - Q*Q'; outer=[N M]; V = outer * inv(middle) * outer';
