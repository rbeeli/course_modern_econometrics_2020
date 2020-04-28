function [pvaluevec, phat]=ButPao(Y,X,c)
% computes the sequential vector of p-values for testing AR(maxp)
%   for the regression model Y=XB+e, where e is presumed AR(p)
% INPUT
% time series column vector Y
% regression matrix X
%  if omitted, default is column of ones
%  Can pass [] for no X matrix.
% c is an maxp-length vector of significance levels, with default
%   maxp=7 and c=[c_1,..,c_maxp]=[0.175 0.15 0.10 0.075 0.05 0.025]
% OUTPUT: 
% pvaluevec is the vector of p-values, starting with AR(1).
% phat is the estimated AR(p) order, based on the p-values, and c.

global Omega G T k r maxp
if nargin<3, maxp=7; c=[0.175 0.15 0.125 0.10 0.075 0.05 0.025]; else maxp=length(c); end
T=length(Y); if nargin<2, X=ones(T,1); end
if isempty(X), k=0; else [~,k]=size(X); end
pvaluevec=NaN(maxp,1); r=NaN(maxp,1);
if isempty(X), M=eye(T); else M=makeM(X); end 
if isempty(X), G=eye(T); else G=makeG(X); end % G is such that M=G'G and I=GG'
e=M*Y; for i=1:maxp, r(i)=(e'*makeA(T,i)*e)/(e'*e); end
pvaluevec(1)=cdfratio(r(1),G*makeA(T,1)*G',eye(T-k),eye(T-k),[],2);
Omega=eye(T-k); % For pure AR(p) case, null is Identity.
% Omegai=G*Psii*G'; Omega=inv(Omegai); % More general varcov
options = optimset('Display','off'); 
phat=-1; m=1; sinit=0; s=fsolve(@spe,sinit,options);
P=makeP(s); [V,D]=eig(P); D=diag(D); tol=1e-6;
if any(D)<tol, disp('P<=0'), D(D<tol)=tol; P=V*diag(D)*V'; end
Pi=inv(P); H=makeH(m,Pi); [V,D]=eig(H); D=diag(D); tol=1e-6;
if any(D)<tol, disp('H<=0'), D(D<tol)=tol; H=V*diag(D)*V'; end
for m=2:maxp
  Pm1=P; Pim1=Pi; Hm1=H;  
  sinit=zeros(m,1); s=fsolve(@spe,sinit,options);
  P=makeP(s); [V,D]=eig(P); D=diag(D); tol=1e-6;
  if any(D)<tol, disp('P<=0'), D(D<tol)=tol; P=V*diag(D)*V'; end
  Pi=inv(P); H=makeH(m,Pi); [V,D]=eig(H); D=diag(D); tol=1e-6;
  if any(D)<tol, disp('H<=0'), D(D<tol)=tol; H=V*diag(D)*V'; end
  sm=s(end); w0=sign(sm)*sqrt( log( det(P)/det(Pm1) ) );
  v0=sm*sqrt(det(H)/det(Hm1)); v0=v0*( tr(Pim1)/tr(Pi) )^(m-1);
  if (isreal(w0) && isreal(v0)) 
    pvaluevec(m)=normcdf(w0) + normpdf(w0)*(1/w0 - 1/v0);
  end
end
if all(isreal(pvaluevec))
  phat=0;  
  for i=1:maxp
    if (pvaluevec(i)<c(i)) || (pvaluevec(i)>(1-c(i))), phat=i; end
  end
end

function A=makeA(T,m) % A = 0.5 * 1( |i-j| = m)
v=zeros(T,1); v(m+1)=1; A=0.5*toeplitz(v,v');

function G=makeG(X) % G is such that M=G'G and I=GG'
k=size(X,2); % could also use k = rank(X).
M=makeM(X); % M=eye(T)-X*inv(X'*X)*X', where X is size TXk
[V,D]=eig(0.5*(M+M')); % V are eigenvectors, D eigenvalues
e=diag(D);
[e,I]=sort(e); % I is a permutation index of the sorting
G=V(:,I(k+1:end)); G=G';

function f=spe(s),  global G k T r
m=length(s); f=zeros(m,1);
for i=1:m
  Pi=inv(makeP(s)); GAG=G*makeA(T,i)*G';
  f(i)=tr(Pi*(GAG-r(i)*eye(T-k))); 
end
%if 1==1, f=sum(f.^2); end

function P = makeP(s),  global Omega G k T r
m=length(s); Sum=zeros(T-k,T-k); 
for i=1:m
  GAG=G*makeA(T,i)*G'; Sum=Sum+s(i)*GAG; 
end
rr=r(1:m); P = Omega + 2*(rr'*s)*eye(T-k) - 2*Sum;

function H = makeH(m,Pinv),  global G T k r
H=zeros(m,m); I=eye(T-k);
for i=1:m,  Ai=G*makeA(T,i)*G';
  for j=1:m,  Aj=G*makeA(T,j)*G'; 
    H(i,j)=2*tr( Pinv*(Ai-r(i)*I)*Pinv*(Aj-r(j)*I) ); 
  end
end

function t=tr(Z), t=sum(diag(Z));


