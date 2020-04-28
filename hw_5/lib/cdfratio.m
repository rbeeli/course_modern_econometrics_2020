function F=cdfratio(rvec,A,B,Sigma,mu,method)
if nargin<6, method=1; end
if nargin<5 || isempty(mu), mu=zeros(length(A),1); end
[V,D]=eig(0.5*(Sigma+Sigma')); W=sqrt(D); Sighalf = V*W*V'; SI=inv(Sighalf);
A=Sighalf*A*Sighalf; B=Sighalf*B*Sighalf; mu=SI*mu;
rl=length(rvec); F=zeros(rl,1); n=length(A);
for rloop=1:rl
    r=rvec(rloop); [P,Lam] = eig((A-r*B)); Lam=real(diag(Lam)); v=P'*mu; nc=v.^2;
    if method==1
        F(rloop)=myimhof(0,Lam,ones(n,1),nc);
    else
        [~,cdfval] = spaweightedsumofchisquare(2,0,Lam,ones(n,1),nc);
        F(rloop)=cdfval;
    end
end