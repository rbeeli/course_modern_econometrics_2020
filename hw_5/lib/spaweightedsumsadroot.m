function [s,report]=spaweightedsumsadroot(s0,x,a,df,q,lower,upper,boundtype)

%disp('just checking')

if 1==2 % just plot K'(s), use lower and upper bounds, or choose values depending on the application
  ss=0:0.001:0.8*upper; for i=1:length(ss), s=ss(i); v = 1./(1-2*s*a); kp(i) = sum(a .* v .* (df + q .* v)); end  
  plot(ss,kp-x), title ('K''(s) - x'), error('hello')
end  

tol1 = 1e-8;  % how close to the edge of the support of X'AX
tol2 = 1e-12; % how close to zero Kpp can be in the Newton step
tol3 = 1e-12; % when to stop the root search
tol4 = 1e-8;  % how close to zero shat may be
MAXIT = 200; s=0; iter=0; derivOK=1; disc=10; report=0;

while ( (iter<=MAXIT) & (derivOK==1) & (disc>=tol3)  )
  iter=iter+1; bot=1-2*s0*a; v = 1./bot;
  kp = sum(a .* v .* (df + q .* v));
  kpp= 2 * sum(a.^2 .* v.^2 .* (df + 2 * q .* v));
  derivOK = (kpp > tol2);
  if derivOK
    s = s0 - (kp-x) / kpp;
    if (s<lower) | (s>upper)
      switch boundtype
        case 1 % lower=-Inf; 
          low=-2; s0=low+rand*(upper-low)*0.999;
        case 2 % upper=Inf
          up=10; s0=lower+0.001+rand*(up-lower);
        case 3
          s0=lower+0.001+rand*(upper-lower)*0.999;
        otherwise
          error('boundtype out of range')
      end
      disc=10; 
      %disp(['value of s out of range, trying start ', num2str(s0)])
    else
      disc = abs(s-s0);
      s0 = s;
    end
  else
    report=6; break
  end
end
if report==0
  if     s<lower+tol1, report = 11;
  elseif s>upper-tol1, report = 12;
  elseif iter==MAXIT, report = 3;
  elseif abs(s) < tol4, report = 5; end
end
s=real(s); 

