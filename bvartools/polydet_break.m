function [dyc, dyt]=polydet_break(Y,ck, tt,fig)

% deterministically breaking detrending  with  polynomials
% regressions:
% y(t)= c+ a1*t +a2*t^2+...+ a4*t^4+ e(t) if  t< tt
% y(t)= c+b a1b*(t-tt) +a2b*(t-tt)^2+...+ a4b*(t-tt)^4+ e(t) if  t>= tt(t)

% tt= break date
% ck= order of  the  polynomial
% dyc = estimated cycle
% dyt = estimated trend

enddT=length(Y);
T = 1:1:tt;
Tb = tt+1:1:enddT;
T2= T.^2;
T2b= Tb.^2;
T3= T.^3;
T3b= Tb.^3;
T4= T.^4;
T4b= Tb.^4;
% regressors
if  ck==1
    X = [ones(tt,1), T'];
    Xb = [ones(enddT-tt,1), Tb'];
elseif ck==2
    X = [ones(tt,1), T', T2'];
    Xb = [ones(enddT-tt,1), Tb' T2b'];
elseif ck==3
    X = [ones(tt,1), T', T2', T3'];
    Xb = [ones(enddT-tt,1), Tb' T2b' T3b'];
elseif ck==4
    X = [ones(tt,1), T', T2', T3', T4'];
    Xb = [ones(enddT-tt,1), Tb' T2b' T3b' T4b'];
end

dyo=zeros(enddT,size(Y,2)); dyoo=zeros(enddT,size(Y,2));

for qq=1:size(Y,2)
     dyo(1:tt,qq) = Y(1:tt,qq) - X*((X'*X)\(X'*Y(1:tt,qq))); %cycle
     dyoo(1:tt,qq) = X*((X'*X)\(X'*Y(1:tt,qq)));        % trend
     dyo(tt+1:enddT,qq) = Y(tt+1:enddT,qq) - ...
      Xb*((Xb'*Xb+0.01*eye(size(Xb'*Xb)))\(Xb'*Y(tt+1:enddT,qq))); % cycle
     dyoo(tt+1:enddT,qq) = Xb*((Xb'*Xb+0.01*eye(size(Xb'*Xb)))\(Xb'*Y(tt+1:enddT,qq)));      % trend
     
     if fig==1
         subplot(2,1,1)
          plot(Y(:,qq),'r-','linewidth',2);hold  on; 
          plot(dyoo(:,qq),'k--','linewidth',2);hold  off; axis tight;
          legend('data', 'Poly break trend')
        subplot(2,1,2)
         plot(dyo(:,qq),'b-','linewidth',2); axis tight;
         legend('Poly break cycle')
        pause
     end
end
dyc=dyo;
dyt=dyoo;

end