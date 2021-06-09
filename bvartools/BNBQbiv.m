function [BNperc, BNtra, BQperc, BQtra]=BNBQbiv(Y,nlags,ck,ccorr,stat,fig,spp,timeplot,varname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beveridge-Nelson and Blancahrd-Quah decomposition, see Canova(2007)

% nlags= numbers  of  lags  in VAR
% ck if =1 constant; if =0  no constant
% corr if=1 shocks correlated;  if=0  shocks  uncorrelated
% stat if =1 second variable I(1); if=0 second variable  stationary
% spp  if=1 plot  spectral densities  of  the components

% BNperc  and  BNtra: permanent and  transitory components BN
% BQperc  and  BQtra: permanent and  transitory components BQ

% Fabio Canova,
% latest revision 28-01-2020: fixed loading  on BQ permanent  and  a  few
%                             bugs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


enddT=length(Y);
ddatax1=squeeze(Y(2:enddT,1)-Y(1:enddT-1,1)); % nonstationary 1 variable
time = 1:1:enddT;
if nargin > 7
    time =timeplot;
end
if nargin < 9
    varname = {'log output'};
end



if  stat==0
    ddatax2=Y(2:enddT,2);                         % stationary 2 variable
else
    ddatax2=squeeze(Y(2:enddT,2)-Y(1:enddT-1,2)); % nonstationary 2 variable
end

ddata=[ddatax1 ddatax2];
ma1=mean(ddatax1);
ma2=mean(ddatax2);
ny = size(ddata,2);

if ck==0
    % [Y1,X1] = VAR_str(ddata-mean(ddata),ck,nlags);  % creates X and the Y
    % variables ordered as:  Y(-1), Y(-2), etc,
    [Y1,X1] = YXB_(ddata-mean(ddata),nlags,[0 0]);
else
    % [Y1,X1] = VAR_str(ddata,ck,nlags);  % creates X and the Y
    [Y1,X0] = YXB_(ddata,nlags,[ck 0]);
    % variables ordered as:  constant, Y(-1), Y(-2), etc.
    X1 = X0(:,[ny*nlags+1,1:ny*nlags]);
end

% creates Y and lags
Bet1=(X1'*X1)\(X1'*Y1);  % beta ols
%Bet1=(X1'*X1+0.001*eye(length(X1'*X1)))\(X1'*Y1);  % beta ridge
res1=Y1-X1*Bet1;      % residuals ols
Sigma1=cov(res1);   % sigma ols

css1=zeros(2,2);
for q=1:nlags
    css1(1,1)=css1(1,1)+Bet1(2*(q-1)+ck+1,1);
    css1(2,1)=css1(2,1)+Bet1(2*(q-1)+ck+2,1);
    css1(1,2)=css1(1,2)+Bet1(2*(q-1)+ck+1,2);
    css1(2,2)=css1(2,2)+Bet1(2*(q-1)+ck+2,2);
end
% Bet1=estimated coefficients; css1= sum of  coefficients
lr1=(eye(2)-css1)\eye(2);   % LR multiplier (equivalent to inv(eye(2)-css1)

if  ck==1
    cb1(1)=Bet1(1,1); cb1(2)=Bet1(1,2); % intercepts
    dd1=((eye(2)-css1)\cb1')';  % constants
end

% restriction  is D(1)A_0 A_0' D(1)= lr *Sigma1*lr'
% so D(1)*A_0 is  cholesky of (lr*Sigma1*lr')
chlr1=chol(lr1*Sigma1*lr1')';  % cholesky of  LR multiplier
Q=(lr1 \ eye(size(lr1)))* chlr1;
res=res1/Q;                    % structural  shocks

% BQ decomposition for correlated  shocks: Cover et  al  (2003)
aalp=-css1(1,2)/(1-css1(2,2));
BQlr1=[1.0/(1+aalp) aalp/(1+aalp)
    -1.0/(1+aalp)  1.0/(1+aalp)]\eye(2);

% permanent component
BQper=zeros(enddT,2); BNper=zeros(enddT,2);
if ccorr==0
    if ck==0
        BNper(1+nlags:enddT-1,1)=ma1(1)+ lr1(1,1)*res1(1:enddT-nlags-1,1)+ ...
            lr1(1,2)*res1(1:enddT-nlags-1,2);
        BQper(1+nlags:enddT-1,1)=ma1(1)+ (lr1(1,1)*Q(1,1))*res(1:enddT-nlags-1,1);  %res(:,1)=supply shock
        if  stat==0
            BNper(1+nlags:enddT-1,2)=ma2(1);
            BQper(1+nlags:enddT-1,2)=ma2(1);
        else
            BNper(1+nlags:enddT-1,2)=ma2(1)+ lr1(2,1)*res1(1:enddT-nlags-1,1)+ ...
                lr1(2,2)*res1(1:enddT-nlags-1,2);
            BQper(1+nlags:enddT-1,2)=ma2(1)+ (lr1(2,1)*Q(2,1))*res(1:enddT-nlags-1,1);  %res(:,1)=supply shock
        end
    else
        BNper(1+nlags:enddT-1,1)=dd1(1)+ lr1(1,1)*res1(1:enddT-nlags-1,1)+ ...
            lr1(1,2)*res1(1:enddT-nlags-1,2);
        BQper(1+nlags:enddT-1,1)=dd1(1)+ (lr1(1,1)*Q(1,1))*res(1:enddT-nlags-1,1);
        if  stat==0
            BNper(1+nlags:enddT-1,2)=dd1(2);
            BQper(1+nlags:enddT-1,2)=dd1(2);
        else
            BNper(1+nlags:enddT-1,2)=dd1(2)+ lr1(2,1)*res1(1:enddT-nlags-1,1)+ ...
                lr1(2,2)*res1(1:enddT-nlags-1,2);
            BQper(1+nlags:enddT-1,2)=dd1(2)+ (lr1(2,1)*Q(2,1))*res(1:enddT-nlags-1,1);
        end
    end
else
    if ck==0
        BNper(1+nlags:enddT-1,1)=ma1(1)+ lr1(1,1)*res1(1:enddT-nlags-1,1)+ ...
            lr1(1,2)*res1(1:enddT-nlags-1,2);
        BQper(1+nlags:enddT-1,1)=ma1(1)+ BQlr1(1,1)*res1(1:enddT-nlags-1,1)+ ...
            BQlr1(1,2)*res1(1:enddT-nlags-1,2);
        if  stat==0
            BNper(1+nlags:enddT-1,2)=ma2(1);
            BQper(1+nlags:enddT-1,2)=ma2(1);
        else
            BNper(1+nlags:enddT-1,2)=ma2(1)+ lr1(2,1)*res1(1:enddT-nlags-1,1)+ ...
                lr1(2,2)*res1(1:enddT-nlags-1,2);
            BQper(1+nlags:enddT-1,2)=ma2(1)+BQlr1(2,1)*res1(1:enddT-nlags-1,1)+ ...
                BQlr1(2,2)*res1(1:enddT-nlags-1,2);
        end
    else
        BNper(1+nlags:enddT-1,1)=dd1(1)+ lr1(1,1)*res1(1:enddT-nlags-1,1)+ ...
            lr1(1,2)*res1(1:enddT-nlags-1,2);
        BQper(1+nlags:enddT-1,1)=dd1(1)+ BQlr1(1,1)*res1(1:enddT-nlags-1,1)+ ...
            BQlr1(1,2)*res1(1:enddT-nlags-1,2);
        if  stat==0
            BNper(1+nlags:enddT-1,2)=dd1(2);
            BQper(1+nlags:enddT-1,2)=dd1(2);
        else
            BNper(1+nlags:enddT-1,2)=dd1(2)+ lr1(2,1)*res1(1:enddT-nlags-1,1)+ ...
                lr1(2,2)*res1(1:enddT-nlags-1,2);
            BQper(1+nlags:enddT-1,2)=dd1(2)+BQlr1(2,1)*res1(1:enddT-nlags-1,1)+ ...
                BQlr1(2,2)*res1(1:enddT-nlags-1,2);
        end
    end
end

% cumulating  permanent components
BQperc=zeros(enddT,2); BNperc=zeros(enddT,2);
BN0=Y(nlags,1); BQ0=Y(nlags,1); % make sure  the  initial conditions match
for q=1+nlags:enddT-1
    BNperc(q,1)= BN0+BNper(q,1);  BN0=BNperc(q,1);  % output
    BQperc(q,1)= BQ0+BQper(q,1);  BQ0=BQperc(q,1);
end
if  stat==0
    BNperc(:,2)= BNper(:,2);    % hours
    BQperc(:,2)= BQper(:,2);
else
    BN0=Y(nlags,2); BQ0=Y(nlags,2); % make sure  the  initial conditions match
    for q=1+nlags:enddT-1
        BNperc(q,2)= BN0+BNper(q,2);  BN0=BNperc(q,2);  % output
        BQperc(q,2)= BQ0+BQper(q,2);  BQ0=BQperc(q,2);
    end
end

% transitory
BQtra=zeros(enddT,2); BNtra=zeros(enddT,2);
for q=1+nlags:enddT-1
    BNtra(q,:)=Y(q,:)-BNperc(q,:);   % level data - level permanent
    BQtra(q,:)=Y(q,:)-BQperc(q,:);
end



if fig==1
    figure(1)
    subplot(2,2,1)
    plot(time(1+nlags:enddT-1),BNperc(1+nlags:enddT-1,1),'k--','linewidth',2); hold  on; ...
        plot(time(1+nlags:enddT-1),Y(1+nlags:enddT-1,1), 'r', 'linewidth',2); hold  off; axis  tight;
    title(['BN ' varname{:}])
    legend('Permanent BN','data')
    subplot(2,2,3)
    plot(time(1+nlags:enddT-1),BNtra(1+nlags:enddT-1,1),'b', 'linewidth',2); axis  tight;
    legend('BN Transitory')
    subplot(2,2,2)
    plot(time(1+nlags:enddT-1),BQperc(1+nlags:enddT-1,1),'k--','linewidth',2); hold  on; ...
        plot(time(1+nlags:enddT-1),Y(1+nlags:enddT-1,1), 'r', 'linewidth',2); hold  off; axis  tight;
    title(['BQ ' varname{:}])
    legend('Permanent BQ','data')
    subplot(2,2,4)
    plot(time(1+nlags:enddT-1),BQtra(1+nlags:enddT-1,1),'b', 'linewidth',2); axis  tight;
    legend('BQ Transitory')
    pause
end

if spp==1
    [ww2, ~]=pwelch(BNtra(:,1));
    [ww3, ~]=pwelch(BQtra(:,1));
    [ww4, ~]=pwelch(Y(:,1));
    [ww5, ~]=pwelch(BNperc(:,1));
    [ww6, ~]=pwelch(BQperc(:,1));
    figure(2)
    subplot(2,1,1)
    plot(log(ww2),'k--', 'linewidth',2); hold on;
    plot(log(ww3),'b-.','linewidth',2); hold  on;
    plot(log(ww4),'r', 'linewidth',2); hold off; axis tight;
    legend('BN transitory', 'BQ transitory', 'data')
    title('log Spectra')
    subplot(2,1,2)
    plot(log(ww5),'k--', 'linewidth',2); hold  on;
    plot(log(ww6),'b-.','linewidth',2); hold on;
    plot(log(ww4),'r', 'linewidth',2); hold off; axis  tight;
    legend('BN permanent', 'BQ permanent', 'data')
    pause
    close all;
end

end

