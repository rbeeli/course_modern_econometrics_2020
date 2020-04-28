function [param, stderr, resid, varcov, loglik]=armareg(y,X,p,q,exact)
% Estimation of Y=XB+e, where e is ARMA(p,q)
% Set exact=1 for exact ML, otherwise conditional ML is used.
% param=[B ; ar terms ; ma terms ; sigma]
% stderr is same shape as param and gives approximate standard errors
% resid is the estimated white noise series
% varcov is the entire (estimated) variance covariance matrix
% Pass X as [] if there is no constant term.
% If X is a scalar, it is set to a vector of ones

ylen=length(y); y=reshape(y,ylen,1);
if length(X)==1, X=ones(ylen,1); end

if isempty(X)
  res=y; beta=[]; nrow=ylen; ncol=0; 
else
  [nrow,ncol]=size(X); beta=inv(X'*X)*X'*y; res=y-X*beta; %#ok<MINV>
end
if p+q==0 % just use OLS and exit
  sigma=sqrt(res'*res/ylen); param=[beta' sigma]'; return
end  

% % use Matlab's routine for ARMA estimation to get initial values. VERY fast!
% m=armax(res,[p,q]); phiinit=-m.a(2:end); thetainit=m.c(2:end);
% if p>0 % if the estimated AR poly is nonstationary, just set the values to zero.
%   rootcheck=min(abs(roots(m.a(end:-1:1)))); if rootcheck<=1, phiinit=zeros(1,p); end
% end
% if q>0 % if noninvertible, just set the values to zero.
%   rootcheck=min(abs(roots(m.c(end:-1:1)))); if rootcheck<=1, thetainit=zeros(1,q); end
% end
% initvec=[beta' phiinit thetainit sqrt(m.noise)]';

initvec=[beta' zeros(1,p+q) std(y)]';

mytol=1e-3; ablo=-1+mytol; abhi=1-mytol; sclo=mytol;
if (p+q)==1 % for an AR(1) or MA(1) model. No restrictions on the regression coefficients
  %%%%%%%%        beta        a or b  scale
  bound.lo=   [-ones(1,ncol)   ablo     sclo     ]';
  bound.hi=   [ ones(1,ncol)   abhi     2*std(y) ]';
  bound.which=[zeros(1,ncol)    1       1        ]';
elseif (p==1) && (q==1)
  %%%%%%%%        beta          a     b     scale
  bound.lo=   [-ones(1,ncol)   ablo   ablo     sclo     ]';
  bound.hi=   [ ones(1,ncol)   abhi   abhi     2*std(y) ]';
  bound.which=[zeros(1,ncol)    1     1        1        ]';
else
  bound.which=zeros(1,length(initvec)); % no bounds at all.
end

mletol=1e-4; MaxIter=400; MaxFunEval=MaxIter*length(initvec);
opt=optimset('Display','None','TolX',mletol, ...
    'MaxIter',MaxIter,'MaxFunEval',MaxFunEval,'LargeScale','off');

[pout,negloglik,~,~,~,hess] = ...
   fminunc(@arma_,einschrk(initvec,bound),opt,y,X,p,q,exact,bound);
loglik=-negloglik; varcov=inv(hess); 
% now transform back and apply delta method to get varcov
[param,varcov]=einschrk(pout,bound,varcov); 

if nargout>1 % get varcov and standard errors
  if 1==2 % direct Hessian calc instead of bfgs output
    H = -hessian(@arma_,param,y,X,p,q,exact); varcov=inv(H);
  end  
  stderr=sqrt(diag(varcov));
end

if nargout>2 % get residuals
  littlesig=param(end);  
  if exact==1
    if isempty(X), z=y; 
    else beta = param(1:ncol); z=y-X*beta; end
    a=param(ncol+1:ncol+p); b=param(ncol+p+1:end-1);
    Sigma = acvf(a,b,nrow); SigInv=inv(Sigma);
    [V,D]=eig(0.5*(SigInv+SigInv')); W=sqrt(D); SigInvhalf = V*W*V';
    resid = SigInvhalf*z/littlesig;
  else
    [~,uvec]=arma_(param,y,X,p,q,0); resid=uvec/littlesig;
  end  
end  



function [loglik,uvec]=arma_(param,y,X,p,q,exact,bound)

if nargin<7, bound=0; end
if isstruct(bound)
  param=einschrk(real(param),bound,999);
end
if any(isinf(param)) || any(isnan(param))
  param=zeros(length(param),1); param(end)=1;
end

if isempty(X), nrow=length(y); ncol=0;
else [nrow,ncol]=size(X); regbeta=param(1:ncol); end  

a=param(ncol+1:ncol+p); b=param(ncol+p+1:end-1); 
sig=abs(param(end)); % this is NOT sigma^2, but just (little) sigma.

if p>0 % enforce stationarity
  theroots = roots([ -a(end:-1:1); 1 ]);
  rootcheck=min(abs(theroots));
  if rootcheck<=1.0001, loglik=abs(1.01-rootcheck)*1e6; return, end
end
if q>0 % enforce invertibility
  rootcheck=min(abs(roots([ b(end:-1:1); 1 ])));
  if rootcheck<=1.0001, loglik=abs(1.01-rootcheck)*1e6; return, end
end

if isempty(X), z=y; else z=y-X*regbeta; end

if exact==1 % get the exact likelihood
  uvec=0;
  if (p==1) && (q==0) % speed this case up considerably.
    K=(-nrow/2) * log(2*pi); 
    s2=sig^2;
    e=z(1)^2*(1-a^2) + sum( (z(2:end) - a*z(1:end-1)).^2 );
    % True likelihood is: 
    % ll = 0.5*log(1-a^2) + K - nrow * log(sig) - e/2/s2; 
    % but if you include K, then this is not compatible with 
    % the general case below, and so LRTs would be wrong.
    %ll =     0.5*log(1-a^2) - nrow * log(sig) - e/2/s2;
    ll = K + 0.5*log(1-a^2) - nrow * log(sig) - e/2/s2;
  else    
    Sigma=acvf(a,b,nrow); Vi=inv(Sigma); detVi=det(Vi);
    if detVi<=0, loglik=abs(detVi+0.01)*1e4; return, end
    ll = -nrow * log(sig) + 0.5*log(detVi) - z'*Vi*z/(2*sig^2); %#ok<MINV>
  end  
else % conditional likelihood
  reversearvec= a(p:-1:1); % avoid reversing the part of z each time
  uroll=zeros(q,1); % a rolling window of U_t hat values
  uvec=zeros(nrow-p,1); % all the T-p U_t hat values
  for t=p+1:nrow
    u=z(t);
    if p>0, u=u-sum( z((t-p):t-1).*reversearvec ); end
    if q>0, u=u-sum(uroll.*b); uroll=[u ; uroll(1:q-1)]; end
    uvec(t-p)=u;
  end
  ll = - nrow * log(sig) - sum(uvec.^2)/(2*sig.^2);
end
loglik = -ll;
