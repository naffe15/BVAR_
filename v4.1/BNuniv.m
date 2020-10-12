function [yc,  yt]=BNuniv(Y,nlags,ck,fig,lr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beveridge-Nelson univariate decomposition, see Canova(2007)

% nlags number  of  lags  in  AR regression
% ck if =1 constant;  if  =0 no  constant
% lr if =1 use  sample mean; if =0  use  long run mean

% yc = estimated cycle
% yt = estimated trend

% Fabio Canova,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

enddT=length(Y);
yc=[]; yt=[];

for  qq=1:size(Y,2)
    ddatax=squeeze(Y(2:enddT,qq)-Y(1:enddT-1,qq));
    %[Z1,W1] = VAR_str(ddatax-mean(ddatax),ck,nlags);  % creates Y and lags
    % variables ordered as:  Y(-1), Y(-2), etc, constant,
    [Z1,W0] = YXB_(ddatax-mean(ddatax),nlags,[ck 0]);  % creates Y and lags
    % variables ordered as:  constant, Y(-1), Y(-2), etc.
    W1 = W0(:,[nlags+1,1:nlags]);
    Bet1x=(W1'*W1)\(W1'*Z1);  % beta ols
    %Bet1x=(W1'*W1+0.001*eye(length(W1'*W1)))\(W1'*Z1);  % beta ridge
    res1x=Z1-W1*Bet1x;      % residuals ols
    if ck==1
        ccx1=Bet1x(1,1);
    else
        ccx1=0.0;
    end
    ssx1=0.0;
    for q=1:nlags
        ssx1=ssx1+Bet1x(q+ck,1);
    end
    lrx1=1.0/(1-ssx1);        % LR multiplier
    
    if lr==1
        ddx1=mean(ddatax(:,1));
    else
        ddx1=lrx1*ccx1;
    end
    
    % permanent component
    BNpotx  = zeros(enddT,1);
    BNpotxc = zeros(enddT,1);
    BN0     = Y(nlags,qq);
    BNpotx(nlags+1:enddT-1,1) = ddx1 + lrx1*res1x(1:enddT-nlags-1,1);
    % cumulating  permanent component
    for q=nlags+1:enddT
        BNpotxc(q,1) = BN0 + BNpotx(q,1);
        BN0          = BNpotxc(q,1);
    end
    % transitory component
    BNgapx=zeros(enddT,1);
%    BNgapx(nlags+1:enddT,1)=Y(nlags+1:enddT,qq)-BNpotxc(1:enddT-nlags,1);
    BNgapx(nlags+1:enddT,1)=Y(nlags+1:enddT,qq)-BNpotxc(nlags+1:enddT,1);
    
    %[Y(:,qq) BNpotxc BNgapx]
    
    if fig==1
        subplot(2,1,1)
        plot(BNpotxc(nlags+1:enddT), 'k--', 'linewidth',2); hold on;
        plot(Y(nlags+1:enddT,qq), 'r', 'linewidth',2); hold off;
        axis  tight, legend('BN univariate trend', 'data');
        subplot(2,1,2)
        plot(BNgapx(nlags+3:enddT-1), 'b', 'linewidth',2); %hold on;
        %plot(Y(:,qq), 'b', 'linewidth',2); hold off;
        axis  tight, legend('BN univariate cycle');
        pause
    end
    
    yc=[yc BNgapx];
    yt=[yt BNpotxc];
end

end