function [dyc, dyt]=hamfilter(X,h,d,ck,fig,timeplot, nameplot)
%HAMFILTER Hamilton (2018) regression filter.
%
%   [dyc,dyt] = HAMFILTER(X,h,d,ck,fig)
%   [dyc,dyt] = HAMFILTER(X,h,d,ck,fig,timeplot,nameplot)
%
%   For each series, regresses y(t) on y(t-h), ..., y(t-h-d+1); the
%   fitted value is the trend and the residual the cycle.
%
%   Inputs
%     X         T x n    series, one per column
%     h         scalar   horizon of the projection (e.g. 8 for quarterly)
%     d         scalar   number of lags (e.g. 4 for quarterly)
%     ck        1 or 0   constant in the regression
%     fig       1 or 0   1 plots each series (pauses after each)
%     timeplot  T x 1    dates for the plots (default 1:T)
%     nameplot  cell     names of the series
%
%   Outputs
%     dyc       T x n    cycle; zeros in the first d+h-1 periods
%     dyt       T x n    trend; zeros in the first d+h-1 periods
%
%   See also POLYDET, BNUNIV, HPFILTER.
%
%   Reference: Hamilton, J. D. (2018), "Why You Should Never Use the
%   Hodrick-Prescott Filter", Review of Economics and Statistics 100(5),
%   831-843.

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
    for v = 1 : size(X,2)
        eval(['titleplot{'   num2str(v) '} =  ''Var' num2str(v) ''';'])
    end
end



for qq=1:QQ
    yh=squeeze(X(d+h:enddT,qq)); %    independent  variable
    R=[];   % reset for each series (old: R kept the lags of the previous series when ck = 0)
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