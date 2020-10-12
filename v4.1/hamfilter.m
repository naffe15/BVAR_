function [dyc, dyt]=hamfilter(X,h,d,ck,fig,timeplot, nameplot)

% local projections (direct  forecast) to  perform detrending
% regressions
% y(t+h)= a*y(t)+b*y(t-1)+c*y(t-2)+...+d*y(t-d)+ e(t+h)
% cycle  is  e(t+h)
% trend  is  hat(a)*y(t)+hat(b)*y(t-1)+hat(c)*y(t-2)+...+hat(d)*y(t-d)

% h= horizon of  the  projection
% d= number  of  lags  used
% ck if=1 constant  if =0 no constant in  the projection

% dyc = estimated cycle
% dyt = estimated trend

[enddT,QQ]=size(X);
yc=zeros(enddT,QQ); yt=zeros(enddT,QQ);
R=[];
T = 1:1:enddT;
time = T;
if nargin > 5
    time =timeplot;
end
if nargin > 6
    titleplot = nameplot;
    if length(nameplot) ~= size(X,2)
        error('nameplot size shold be the same as the column of Y')
    end
else
    for v = 1 : size(Y,2)
        eval(['titleplot{'   num2str(v) '} =  ''Var' num2str(v) ''';'])
    end
end



for qq=1:QQ
    yh=squeeze(X(d+h:enddT,qq)); %    independent  variable
    if  ck==1
        R=ones(enddT-d-h+1,1);     % constant
    end
    
    for  jj=1:d
        r=squeeze(X(d+1-jj:enddT-h+1-jj,qq));
        R=[R r];  % dependent  variables
    end
    
    yc(d+h:enddT,qq) = yh - R*((R'*R)\(R'*yh)); % cycle
    yt(d+h:enddT,qq) = R*((R'*R)\(R'*yh));        % trend
    
    if  fig==1
        subplot(2,1,1)
        %plot(time(d+h:enddT),yh(1:enddT-d-h,1),'r', 'linewidth',2);hold  on;
        plot(time(d+h:enddT),X(d+h:enddT,qq),'r', 'linewidth',2);hold  on;
        plot(time(d+h:enddT),yt(d+h:enddT,qq),'k--','linewidth',2);hold off; axis  tight;
        legend('data', 'Hamil trend')
        title(titleplot(qq))
        subplot(2,1,2)
        plot(time(d+h:enddT),yc(d+h:enddT,qq),'b', 'linewidth',2);
        legend('Hamil cycle')
        pause
    end
end
dyc=yc;
dyt=yt;
end