function cdf = myimhof (xvec,wgt,df,nc,uplim)
% cdf = myimhof (xvec,wgt,df,nc,uplim)
% default is to use the transformation to (0,1). If you do not 
%   want this, pass uplim>0, so integration is from (0,uplim), without the
%   transformation.

if nargin<5, uplim=-1; end

k=length(wgt); tol=1e-9;
if nargin < 3, df=ones(k,1); end
if nargin < 4, nc=zeros(k,1); end

lx=length(xvec); cdf=zeros(lx,1);
for i=1:lx
  x=xvec(i);
  if 1==1 % Use the quadgk routine in matlab; see its help file. 
      % it allows for singularities at the end points, so do the transformation
    cdf(i) = 0.5 - (1/pi) * quadgk(@(tvec)ff(tvec,x,wgt,df,nc,-99),0,1);
  else % use quadl
    if uplim<0, cdf(i) = 0.5 - (1/pi) * quadl(@ff,tol,1-tol,tol,0,x,wgt,df,nc,uplim);
    else cdf(i) = 0.5 - (1/pi) * quadl(@ff,tol,uplim,tol,0,x,wgt,df,nc,uplim);
    end
   
  end
end  

function I = ff(tvec,x,wgt,df,nc,uplim)
vlen=length(tvec); I=zeros(1,vlen);
usetransform=(uplim<0);
for ii=1:vlen
  t=tvec(ii); 
  if usetransform, u=(1-t)/t; else  u=t; end
  p = u*wgt; b=p.^2; c=1+b;
  ss0=sum(df.*atan(p) + nc.*p./c); beta = (ss0 - x*u)/2;
  ss1=sum(nc.*b./c); ss2=sum(df.*log(c)); gam=exp(0.5*ss1 + 0.25*ss2);
  I(ii) = sin(beta) / (u*gam);
  if usetransform, I(ii) = I(ii) / t^2; end % don't forget du/dt
end
