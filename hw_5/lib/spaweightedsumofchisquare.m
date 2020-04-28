function [pdf,cdf,svec,K,Kpp,kap3,forboth] = spaweightedsumofchisquare(daniels,xvec,a,df,q,s)
% [pdf,cdf,svec,K,Kpp,kap3,forboth] = spaweightedsumofchisquare(daniels,xvec,a,df,q,s);
% SPA of sum a_i X_i, X_i ~ chi2(df_i,q_i) at elements of x, i.e., q_i is noncentrality param
% pass s as a starting value for shat. 
% set daniels to 2 to use 2nd order approx.
% If you pass s as a scalar, then it is used as a starting value just for xvec(1),
%   and the resulting s^hat is used for xvec(2), and so on.
% Or pass s as a vector the same length as xvec. (THIS IS NOT IMPLEMENTED YET!)
%
% K,Kpp,kap3,forboth correspond to the last element in xvec, i.e., is only useful 
% when calling with a single x value, as I do from sparatio.m

%persistent sstart
%if isempty(sstart), sstart=0; end

if nargin<6, s=0; end
if nargin<5, q = zeros(length(a),1); end
if nargin<4, df = ones(length(a),1); end

sstart=s;

kap3=0; forboth=0; % these are assiged below if daniels=2, otherwise not, but matlab complains if they are not assigned.

alo=2*min(a) ; ahi=2*max(a);
if alo>0, lower=-Inf; upper = 1/ahi; boundtype=1;
  elseif ahi<0, lower=1/alo; upper=Inf; boundtype=2;
  else lower=1/alo; upper=1/ahi; boundtype=3;
end 

pdf=zeros(length(xvec),1); cdf=pdf; svec=pdf;
for i=1:length(xvec)
  x=xvec(i);
  [s,report]=spaweightedsumsadroot(sstart,x,a,df,q,lower,upper,boundtype);
  if (report>0) && (boundtype==3) % try once again, using a grid of s values
    ss=lower:(upper-lower)/200:upper; 
    %x
    for j=2:length(ss)-1
      [s,report]=spaweightedsumsadroot(ss(j),x,a,df,q,lower,upper,boundtype); 
      if report==0
        %disp('valid saddlepoint found')
        %disp(s)
        break
      end
    end
  end  
  svec(i)=s;
  sstart=s; % turn on to use shat from xvec(i-1) as starting value
  
  if ((report==0) || (report==5))
    v = 1./(1-2*s*a);
    K = 0.5 * sum(df .* log(v)) + s*sum(a .* q .* v);
    Kpp  = 2 * sum(a.^2 .* v.^2 .* (df + 2 * q .* v));
    pdf(i) = exp(K - x*s) / sqrt(2*pi*Kpp);
    fac = 2*s*x-2*K; ww=sign(s)*sqrt(fac); u=s*sqrt(Kpp);
    if ((report==5) || (abs(ww)<1e-7))
      %disp(['x = ',num2str(x),' is (or is very close to) the expected value'])
      Kpp0  = 2 * sum(a.^2 .* (df + 2 * q));
      Kppp0 = 8 * sum(a.^3 .* (df + 3 * q));
      cdf(i)=0.5+Kppp0/6/sqrt(2*pi)/Kpp0^(3/2);

      K3= 8*sum(a.^3 .* v.^3 .* (df + 3*q.*v));  
      K4=48*sum(a.^4 .* v.^4 .* (df+4*q.*v));
      kap3= K3/(Kpp)^(3/2);  kap4=K4/(Kpp)^(4/2);
      forboth = (kap4/8 - 5*kap3^2/24);
    else
      if (~isreal(ww)) || (~isreal(u))
       cdf(i)=-1; % flag that it failed
       continue
      end
      npdf=normpdf(ww);
      tempcdf=normcdf(ww)+npdf*(1/ww - 1/u);
      if daniels==2
        K3= 8*sum(a.^3 .* v.^3 .* (df + 3*q.*v));  
        K4=48*sum(a.^4 .* v.^4 .* (df+4*q.*v));
        kap3= K3/(Kpp)^(3/2);  kap4=K4/(Kpp)^(4/2);
        forboth = (kap4/8 - 5*kap3^2/24);

        % the cdf correction can be problematic!
        bigterm = forboth/u - 1/u^3 - kap3/2/u^2 + 1/ww^3;
        correction=npdf * bigterm;
        if abs(correction)/tempcdf > 0.1
          % 'The 2nd order CDF correction term is being ditched!'
          correction=0;
        end  
        cdf(i) = tempcdf - correction;
        if ((cdf(i)<0) || (cdf(i)>1))
          %'The 2nd order CDF correction term is being ditched!'
          cdf(i) = tempcdf;
        end  
      else
        cdf(i) = tempcdf;
      end
    end
    if daniels==2, pdf(i) = pdf(i) * (1+forboth); end
  else
    pdf(i)=0; cdf(i)=0;
  end
end

  
%  just for reference:
% v = 1./(1-2*s*a);
% K    = 0.5 * sum(df .* ln(v)) + s*sum(a .* q .* v);
% Kp   =     sum(a    .* v    .* (df +     q .* v));
% Kpp  = 2 * sum(a.^2 .* v.^2 .* (df + 2 * q .* v));
% Kppp = 8 * sum(a.^3 .* v.^3 .* (df + 3 * q .* v));
