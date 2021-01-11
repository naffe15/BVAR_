function [dyc, dyt]=polydet(Y,ck, fig, timeplot, nameplot)

% deterministically  detrending  with  polynomials
% regression: y(t)= c+ a1*t +a2*t^2+...+ a4*t^4+ e(t)
% the cycle  is  e(t)
% the trend is hat(c}+ hat(a1)*t+hat(a2)*t^2 +..+hat(a4)*t^4

% ck= order of  the  polynomial
% dyc = estimated cycle
% dyt = estimated trend

enddT=length(Y);
T = 1:1:enddT;
time = T;
if nargin > 3
    time =timeplot;
end
if nargin >4
    titleplot = nameplot;
    if length(nameplot) ~= size(Y,2)
        error('nameplot size shold be the same as the column of Y')
    end
else
    for v = 1 : size(Y,2)
        eval(['titleplot{'   num2str(v) '} =  ''Var' num2str(v) ''';'])
    end
end


T2= T.^2;
T3= T.^3;
T4= T.^4;
% regressors
if  ck==1
    X = [ones(enddT,1), T'];
elseif ck==2
    X = [ones(enddT,1), T', T2'];
elseif ck==3
    X = [ones(enddT,1), T', T2', T3'];
elseif ck==4
    X = [ones(enddT,1), T', T2', T3', T4'];
end

dyo=zeros(enddT,size(Y,2)); dyoo=zeros(enddT,size(Y,2));

for qq=1:size(Y,2)
    dyo(:,qq) = Y(:,qq) - X*((X'*X)\(X'*Y(:,qq))); % cycle
    dyoo(:,qq) = X*((X'*X)\(X'*Y(:,qq)));          % trend
    
    if fig==1
        subplot(2,1,1)
        plot(time,Y(:,qq),'r-','linewidth',2);hold  on;
        plot(time,dyoo(:,qq),'k--','linewidth',2);hold  off; axis tight;
        legend('data', 'Polynomial trend')
        title(titleplot(qq))
        subplot(2,1,2)
        plot(time,dyo(:,qq),'b-','linewidth',2); axis tight;
        legend('Polynomial cycle')
        pause
    end
end
dyc=dyo;
dyt=dyoo;

end