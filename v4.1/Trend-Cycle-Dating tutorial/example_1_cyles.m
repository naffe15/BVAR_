clear all
close all
clc

addpath ../../cmintools/
addpath ../../v4.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% realfin_exa.m, 21-11-2019, version 3.2, fabio canova
%
% Freely available  for  distribution. Use  at  your  own  risk
%
% this program illustrates the use of different filters to compute
% trend/cycle  decompositions using  euro area data  on real (GDP) and
% financial (credit-to-GDP) variables
%
% data for  credit  to  GDP  available  only  from  1999:1
% so  long  cycles  (> 40 quarters) are  poorly estimated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% exercises  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  exercise 1: Vary the upper and lower limit of BK and CF filters
%             for example use 8 and 24 or  8 to  64 and see what happens.
%  exercise 2: Change the assumed properties of the time series in CF
%             filter and compare results. In particular, check what happens
%             if  the  data is  assumed  to  have  a  unit  root or  not.
%  exercise 3: Change  the  long  differencing  filter  to  8  years
%  exercise 4: Change h from 8 to 12 and d from 4 to  6 in Hamilton filter.
%  exercise 5: Use  shorter  lag  length  in  the  univariate BN
%  exercise 6: Allow  trend  and  cycle  to  be  correlated in  UC
%  exercise 7: Vary parameters  of  the  butterworth filters
%  exercise 8: use  ff=0 or  ka=0 in  BN -BQ  bivaraite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=======================================================================
% LOAD DATA
%=======================================================================
% load the THE AWM DATABASE: Quarterly
[a,b,~] = xlsread('awm19up18.csv');
% names of variables
varnames = b(1,2:end);

% time convention: Q1 = .00 and Q4 =0.75
time = 1970 : .25 : 2017.75;
time_s   = find(time==1999.50);
time_e   = find(time==2017.75);
time_b   = find(time==2007.75);

[e,f,~]=xlsread('ECB_Credit_gaps.xlsx');
varnames2= f(1,2:end);

% collect real GDP, real consumption, real  investment, GDP defl, HICP,
% short and long term interest rate, commodity  prices,
% labor  productivity, total  employment, credit to NFC to  GDP

% names: YER PCR ITR YED HICP STN LTN COMPR LPROD LNN credit

yy = [a(time_s-2:time_e,1) a(time_s-2:time_e,2) a(time_s-2:time_e,4) ...
    a(time_s-2:time_e,7) a(time_s-2:time_e,19) a(time_s-2:time_e,33) ...
    a(time_s-2:time_e,34) a(time_s-2:time_e,35) ...
    a(time_s-2:time_e,42) a(time_s-2:time_e,29) e(189:264,25)];

timeplot = time(time_s-2:time_e);
varnames = {'log GDP','log credit to GDP'};

% data transformations:
% 1. the log of output, consumption,investment
ddata(:,1:3) = log(yy(1 : length(yy),1:3));
% 2. the log/log difference of GDP deflator index and CPI
ddata(:,4:5) = log(yy(1  : length(yy),4:5));
%ddata(:,4:5) = diff(log(yy(2 : length(yy),4:5)))*400;
% 3. the level of short and long term interest rate
ddata(:,6:7) = yy(1 : length(yy)   ,6:7);
% 4. the log/log difference of commodity prices
ddata(:,8) = log(yy(1  :  length(yy)   ,8));
%ddata(:,8) = diff(log(yy(2 :  length(yy) ,8)))*400;
% 5. the taking the log of labor  productivity and  employment
ddata(:,9:10) = log(yy(1 :  length(yy)  ,9:10));
%ddata(:,9) = diff(log(yy(2: length(yy)   ,9)));
%ddata(:,10) = diff(log(yy(2: length(yy)  ,10));
%  the  log of  credit  to  GDP
ddata(:,11) = log(yy(1: length(yy) ,11));

endd=length(ddata);

% pick log output, log credit/GDP
ddd=ddata(:,[1 11]);

disp('----------------------------------------------------')
disp(' 1) polynomial filter')
disp('---------------------------------------------------')
% parameters
ord=3; % order of the polynomial (up to  4)
gg=1;  % if 1 plot  data, trend and cycle
% dy1  has  the  cycles , dt1 the  trends
% [dy1, dt1]=polydet(ddd,ord,gg);
[dy1, dt1]=polydet(ddd,ord,gg,timeplot,varnames);
close all;


disp(' ----------------------------------------------------')
disp(' 2) HP filter')
disp(' ----------------------------------------------------')
% three  options : lambda=1600, lambda=4, lambda=51200
llam1 = 1600;   % standard  cycles
llam2 = 4;      % short cycles
llam3 = 51200;  % long cycles
% dy2-dy4  have  the cyclical components, dt2-dt4 the trends
dy2=ddd(:,1:size(ddd,2))-Hpfilter(ddd(:,1:size(ddd,2)),llam1);
dy3=ddd(:,1:size(ddd,2))-Hpfilter(ddd(:,1:size(ddd,2)),llam2);
dy4=ddd(:,1:size(ddd,2))-Hpfilter(ddd(:,1:size(ddd,2)),llam3);
dt2=Hpfilter(ddd(:,1:size(ddd,2)),llam1);
dt3=Hpfilter(ddd(:,1:size(ddd,2)),llam2);
dt4=Hpfilter(ddd(:,1:size(ddd,2)),llam3);

%plot(dy2(:,2:2),'r', 'linewidth',2); hold on;
%plot(dy3(:,2:2),'b', 'linewidth',2); hold  on;
%plot(dy4(:,2:2),'k', 'linewidth',2); hold  off;
%legend('Hp1600','Hp4','Hp51200')

% one-sided HP filter
dy5=zeros(endd,size(ddd,2)); dt5=dy5;
for qq=1:size(ddd,2)
    % one-sided_hpfilter_serial works also on ddd directly
    % (no need to loop, here  used only  for plotting)
    disc=0; %(number  of discarded observations, usually 0)
    y=ddd(:,qq);
    % y=ddd(:,1:size(ddd,2));
    [dtx,dyx]=one_sided_hpfilter_serial(y,llam1,disc);
    subplot(2,1,1)
    plot(timeplot,y,'r-','linewidth',2); hold  on;
    plot(timeplot,dt2(:,qq),'b--','linewidth',2); hold  on;
    plot(timeplot,dtx,'k-.','linewidth',2); hold  off; axis tight;
    legend('data','2-sided HP','1-sided HP')
    if  qq==1
        title('log GDP')
    else
        title('log credit to  GDP')
    end
    subplot(2,1,2)
    plot(timeplot,dy2(1:endd,qq),'b--','linewidth',2); hold  on;
    plot(timeplot,dyx(1:endd,1),'k-.','linewidth',2); hold  off; axis tight;
    legend('2-sided HP','1-sided HP')
    pause
    dy5(:,qq)= dyx;
    dt5(:,qq)=dtx;
end
close all;

disp('------------------------------------------------------------')
disp(' 3) differencing')
disp(' -----------------------------------------------------------')
% dy6-dy8 have  the cyclical componentsl dt6-dt8 have  the  trends
% dy6  loose  one  observation; dy7 4  observations, dy8 20 observations
dy6=zeros(endd,size(ddd,2));dy7=zeros(endd,size(ddd,2));
dy8=zeros(endd,size(ddd,2));

% one  quarter difference
dy6(2:endd,1:size(ddd,2))=ddd(2:endd,1:size(ddd,2))-ddd(1:endd-1,1:size(ddd,2));
dt6(2:endd,1:size(ddd,2))=ddd(1:endd-1,1:size(ddd,2));

% one year  difference
dy7(5:endd,1:size(ddd,2))=ddd(5:endd,1:size(ddd,2))-ddd(1:endd-4,1:size(ddd,2));
dt7(5:endd,1:size(ddd,2))=ddd(1:endd-4,1:size(ddd,2));

%  long (5  years) differences
dy8(21:endd,1:size(ddd,2))=ddd(21:endd,1:size(ddd,2))-ddd(1:endd-20,1:size(ddd,2));
dt8(21:endd,1:size(ddd,2))=ddd(1:endd-20,1:size(ddd,2));

for qq=1:size(ddd,2)
    subplot(2,1,1)
    plot(timeplot(1:endd),ddd(1:endd,qq),'g', 'linewidth',2); hold on;
    plot(timeplot(1:endd-1),dt6(1:endd-1,qq),'b-', 'linewidth',2); hold on;
    plot(timeplot(1:endd-4),dt7(1:endd-4,qq),'k-.','linewidth',2); hold on;
    plot(timeplot(1:endd-20),dt8(1:endd-20,qq),'r:','linewidth',2);hold  off;
    legend('data','1odtrend','4odtrend','20odtrend')
    if  qq==1
        title('log GDP')
    else
        title('log credit  to GDP')
    end
    subplot(2,1,2)
    plot(timeplot,dy6(1:endd,qq),'b-', 'linewidth',2); hold on;
    plot(timeplot,dy7(1:endd,qq),'k-.','linewidth',2); hold on;
    plot(timeplot,dy8(1:endd,qq),'r:','linewidth',2);hold  off;
    legend('1odcycle','4odcycle','20odcycle')
    pause
end
close all;


disp('-------------------------------------------------------')
disp(' 4) BP filter ')
disp(' ------------------------------------------------------')
% dy9-dy10 have  the cyclical components; dt8-dt9 the  trends
bp1  = 8;  bp2  = 32;   % band pass filter parameters
% Baxter and King
dy9=zeros(endd,size(ddd,2)); dy10=dy9;
for qq=1:size(ddd,2)
%     dycc = bkfilter(ddd(:,qq),bp1,bp2);
    dycc = bkfilter(ddd(:,qq),bp1,bp2);
    dy9(:,qq)=dycc;
end
dt9=ddd-dy9;

% Christiano-Fitzgerald
% the symmetric filter is the same as bk: cfcpi = cffilter(x,bp1,bp2);
% details about the  parameters of  the asymmetric filter see cffilter.m
for qq=1:size(ddd,2)
    dycf = cffilter(ddd(:,qq),bp1,bp2,1,1,0);
    dy10(:,qq)=dycf;
end
dt10=ddd-dy10;

disp(' ----------------------------------------------------------')
disp(' 5) wavelet  filter')
disp(' ----------------------------------------------------------')
% dy11  has  the cycle; dy12  has  the  cycle+low; dt11-dt12 the  trends
% loose 32  observations
gg=0; % if  1,  plot  business cycle ( business cycle + low frequency)

% need  to  apply  filter  to  inflation (not  log price)
%dddd=[ddd(2:length(ddd),1) diff(ddd(:,2)) ddd(2:length(ddd),3:4)];
[dy11, dt11 ,dy12, dt12]=wavefilter(ddd,gg);

for qq=1:size(ddd,2)
    subplot(2,1,1)
    % start  plotting  trends only from  t=32
    plot(timeplot(32:end),ddd(32:endd,qq)/ddd(32,qq),'r-', 'linewidth',2); hold on;
    plot(timeplot(32:end),dt9(32:endd,qq)/dt9(32,qq),'c--', 'linewidth',2); hold on;
    plot(timeplot(32:end),dt10(32:endd,qq)/dt10(32,qq),'b:', 'linewidth',2); hold on;
    plot(timeplot(32:end),dt11(32:endd,qq)/dt11(32,qq), 'k:', 'linewidth',2);hold on;
    axis tight;
    plot(timeplot(32:end),dt12(32:endd,qq)/dt12(32,qq),'g-.','linewidth',2); hold  off;
    legend('data','BK','CF','WaveBC', 'WavelowBC')
    if  qq==1
        title('log GDP')
    else
        title('log  credit  to  GDP')
    end
    subplot(2,1,2)
    plot(timeplot,dy9(1:endd,qq),'c-', 'linewidth',2); hold on; plot(timeplot,dy10(1:endd,qq),'b', 'linewidth',2); hold on; axis tight;
    plot(timeplot,dy11(1:endd,qq), 'k-', 'linewidth',2);hold on;  plot(timeplot,dy12(1:endd,qq),'g-','linewidth',2); hold  off;
%     legend('BK','CF','WaveBC', 'WavelowBC')
    pause
end
close all;

disp(' -----------------------------------------------')
disp(' 6) Hamilton filter')
disp(' -----------------------------------------------')
% model  y(t+h)= a*y(t)+b*y(t-1)+...+ q*y(t-q)+e(t+h)
% e(t+h) is  the  cyclical
% loose  h+d  observations
% parameters
h=8; % number  of  horizons  forward
d=4; %  number  of  lags
gg=1;  %  plot  cyclical
ff=1;  % include  constant in  the regression
[dy13, dt13]=hamfilter(ddd,h,d,ff,gg,timeplot,varnames);
close all;


disp(' ------------------------------------------------------')
disp(' 7) univariate UC  filter')
disp(' ------------------------------------------------------')
% model  y(t)= yT(t)+yc(T)
%        yT(t)=mu+YT(t-1)+ e(t)
%        yC(t)= phi_1*yC(t-1)+phi_2*yC(t-2)+ u(t)
%        e(t)  ~  N(0, sigma^2_tau         rho*sigma_ tau*sigm_c
%        u(t)       0, rho*sigma_tau*sigma_c   sigma^2_c        )
%  -----------------------------------------------------
%  Careful for  some  series rho neq 0 creates  problems
%  ------------------------------------------------------
% loose  2  observations
% parameters
options.figg     = 1;        % 1: plot figures  of  series+trend/cyclical
options.noprint  = 0;        % 0: print  estimates
options.nsims    = 5000;     % MCMC parameter: number  of simulations
options.burnin   = 30000;    % MCMC parameter: number  of  burn-in
options.rhoyes   = 0;        % if =1;  allow  correlation between trend and  cycle

% prior first moments
% prior mean trend drift 
options.mu0 = 0.1;
% prior mean cycle AR1 and AR2 
options.phi0 = [0.2 0.2]';
% prior mean VARIANCE cycle
options.sigc2 = 0.5;
% prior mean VARIANCE trend
options.sigtau2 = 1.5;
% time for the plot
options.time = timeplot;

dy14=zeros(endd,size(ddd,2)); dt14=zeros(endd,size(ddd,2));
for qq=1:size(ddd,2)
    %     ycycle=[];
    %     ytrend=[];
    y               = squeeze(ddd(:,qq));
%      [ytrend,ycycle,out] = uc_(y,2);
    options.varname =varnames(qq);
    [ytrend,ycycle] = uc1_(y,options);
    dy14(:,qq)      = ycycle;
    dt14(:,qq)      = ytrend;
end
close all;


disp(' ------------------------------------------------------')
disp(' 7b) univariate UC I(2) filter')
disp(' ------------------------------------------------------')
% model  
% from the UC of the form 
% a(t)  = ph1 a(t-1) + ... + phip a(t-p) + ea(t ) [transition 1]
% b(t)  = c(t-1)    + b(t-1) + eb(t);             [transition 2]
% c(t)  = c(t-1)             + ec(t);             [transition 3]
% y(t)  = a(t)      + b(t);                       [transition 4]

opts.time = timeplot;
lags  = 2;
% AR1[Cycle] AR2[Cycle] STD(a)[Cycle]  STD(b)[Trend]  STD(c)[Trend] 
opts.ub = [0.6 0.6 3 3 3];
opts.lb = [0.01 0.01 0.05 0.05 0.05];
opts.max_compute = 2;

% dy14=zeros(endd,size(ddd,2)); dt14=zeros(endd,size(ddd,2));
for qq=1:size(ddd,2)
    %     ycycle=[];
    %     ytrend=[];
    y               = squeeze(ddd(:,qq));
    if qq==1
        opts.varname = 'log GDP';
    elseif qq==2
        opts.varname = 'log Credit to GDP';
    end
    [ytrend,ycycle,out] = uc2_(y,lags,opts);
%     dy14(:,qq)      = ycycle;
%     dt14(:,qq)      = ytrend;
end
close all;

disp('-------------------------------------------------')
disp(' 8) univariate BN filter')
disp('-------------------------------------------------')
% loose  1+nlags  observations
% parameters
lags= 2;  % number  of  lags  in  the  autoregression (keep  it  small
%  otherwise  estimated  long run  mean  is  wrong).
ff=1;     % include  constant in  the  autoregression
gg=1;     % plot cyclical
mm=1;     % =1 use estimate  mean;  =0 use  long run mean
[dy15, dt15] = BNuniv(ddd,lags,ff, gg, mm);
close all;

disp('------------------------------------------------')
disp(' 9) Butterworth filters')
disp('------------------------------------------------')

disp('low  pass  filter')
n1=4;   % order  of  the  polynomial
cutoff1=0.07;   % cutoff frequency;  above 0.05;
[a1,b1]=butter(n1,cutoff1,'low');
%h = fvtool(b1,a1);  % visualization of  the  magnitude/gain  response
dy16=filtfilt(a1,b1,ddd);  % ARMA filter  with weghts given by  a  and  b

disp('band  pass  filter')
n2=4;   % order  of  the  polynomial
cutoff2=[0.07,0.30];   % cutoff frequency
[a2,b2]=butter(n2,cutoff2);
%h = fvtool(b2,a2);  % visualization of  the  magnitude/gain  response
dy17=filtfilt(a2,b2,ddd);  % ARMA filter  with weghts given by  a  and  b


disp(' high  pass filter')
n3=1;   % order  of  the  polynomial
cutoff3=0.07;   % cutoff frequency  % above  0.05  it  becomes  as BP
[a3,b3]=butter(n3,cutoff3,'high');
%h = fvtool(b3,a3);  % visualization of  the  magnitude/gain  response
dy18=filtfilt(a3,b3,ddd);  % ARMA filter  with weghts given by  a  and  b


for qq=1:size(ddd,2)
    figure
    subplot(2,1,1)
    plot(timeplot,dy16(1:endd,qq),'k--', 'linewidth',2); hold on;
    plot(timeplot,ddd(1:endd,qq),'r-', 'linewidth',2); hold off;
    legend('LP','data')
    if  qq==1
        title('log GDP')
    else
        title('log  credit  to  GDP')
    end
    subplot(2,1,2)
    plot(timeplot,dy17(1:endd,qq),'b', 'linewidth',2); hold on;
    plot(timeplot,dy18(1:endd,qq), 'g-.', 'linewidth',2);hold  off; axis  tight;
    legend('BP','HP')
    
    pause
    close all
end


disp('------------------------------------------------')
disp(' 10) bivariate  BN and  BQ filters')
disp('------------------------------------------------')
% loose  1+nlags  observations
% parameters
% nlags=6; % number  of  lags  in  var
% ff=1;  % if 1 constant  in  var
% ka=1;  % if 1 use the mean  rather than estimated  long run  mean
% % in computing  permanent  components
% gg=1;  % if 1 plot  trend  and  cycles of  output
% mm=0;  % if 1  plot  spectral densities
% rhoyes=0;  % if 1 permanent  and  transitory  shocks  correlated
% ssts  =0;  % if 1 second  variable  also  integrated
% % ok for  BN, not  ok  for  BQ
% % dy19-dt19  BN; dy20-dt20 BQ
% [dt19, dy19, dt20, dy20]=BNBQbiv(ddd,nlags,ff,ka,gg,mm, rhoyes,ssts);

nlags   = 6;     % number of lags in var: needs 6/7 to have something reasonable
ff      = 0;     % if 1 constant  in  var (and do  not  demean  the  data)
rhoyes  = 0;     % if 1 shocks are correlated
ssts    = 1;     % if 1 second  variable  also  integrated
gg      = 1;     % if 1 plot  permanent  and  transitory of  log output
mm      = 0;     % if 1 plot  spectral densities

[dt19, dy19, dt20, dy20]=BNBQbiv(ddd,nlags,ff,rhoyes, ssts, gg,mm,timeplot,{'log output'});


close all;

% collect some  cyclical components of  output and plot them
% LT, HP, HP1s, FOD, BP, Hamil, UC, BN, BW

dates=1999.00:0.25:2017.75;

xxc=[dy1(3:endd,1:size(ddd,2)) dy2(3:endd,1:size(ddd,2))  dy5(3:endd,1:size(ddd,2)) dy6(2:endd-1,1:size(ddd,2)) ...
    dy10(3:endd,1:size(ddd,2)) dy13(3:endd,1:size(ddd,2)) dy14(3:endd,1:size(ddd,2)) dy15(3:endd,1:size(ddd,2)) ...
    dy17(3:endd,1:size(ddd,2))];

figure(1)
plot(dates,dy1(1:endd,1),'r', dates,dy2(1:endd,1),'b',dates, dy6(1:endd,1),'g', dates, dy10(1:endd,1),'k', ...
    dates, dy13(1:endd,1),'m', dates, dy14(1:endd,1),'c',dates, dy15(1:endd,1),'y', 'linewidth',2) ; axis  tight;
legend('LT', 'HP', 'FOD', 'BP', 'HAM', 'UC', 'BN')
title('Summary: log output')
pause


figure(2)
plot(dates,dy1(1:endd,2),'r', dates,dy2(1:endd,2),'b',dates, dy6(1:endd,2),'g', dates, dy10(1:endd,2),'k', ...
    dates, dy13(1:endd,2),'m', dates, dy14(1:endd,2),'c',dates, dy15(1:endd,2),'y', 'linewidth',2); axis tight;
legend('LT', 'HP', 'FOD', 'BP', 'HAM', 'UC', 'BN')
title('Summary: log credit to  GDP')
pause

close all;

% autocorrelations and cross correlations of cyclical GDP
% (the dimension of auto1 is 7*(4+1) times 7 with blocks of
% cross correlations at 0,1,2,3,4 (each block is (5 x5)).

xxo1=[dy1(3:endd,1) dy2(3:endd,1) dy6(2:endd-1,1) dy10(3:endd,1) dy13(3:endd,1) dy14(3:endd,1) dy15(3:endd,1) dy17(3:endd,1)];
auto1=autocor(xxo1,4);

disp('contemporaneous correlations of cyclical output')
auto1(1:8,1:8)
disp('first lagged correlations of  cyclical  output')
auto1(9:16,1:8)

% calculation of variabilities
xxg1=[var(dy1(3:endd,1)) var(dy2(3:endd,1)) var(dy6(2:endd-1,1)) ...
    var(dy10(3:endd,1))  var(dy13(3:endd,1)) var(dy14(3:endd,1)) ...
    var(dy15(3:endd,1)) var(dy17(3:endd,1))];
disp('variabilities of  cyclical output')
xxg1(1,:)

% autocorrelations and cross correlations of cyclical credit  to GDP
% (the dimension of auto1 is 7*(4+1) times 7 with blocks of
% cross correlations at 0,1,2,3,4 (each block is (5 x5)).

xxo2=[dy1(3:endd,2) dy2(3:endd,2) dy6(2:endd-1,2) dy10(3:endd,2) ...
    dy13(3:endd,2) dy14(3:endd,2) dy15(3:endd,2) dy17(3:endd,2)];
auto2=autocor(xxo2,4);

disp('contemporaneous correlations of cyclical credit  to  output')
auto2(1:8,1:8)
disp('first lagged correlations of  cyclical credit  to output')
auto2(9:16,1:8)

% calculation of variabilities
xxg2=[var(dy1(3:endd,2)) var(dy2(3:endd,2)) var(dy6(2:endd-1,2)) ...
    var(dy10(3:endd,2))  var(dy13(3:endd,2)) var(dy14(3:endd,2)) ...
    var(dy15(3:endd,2)) var(dy17(3:endd,2))];
disp('variabilities of  cyclical credit  to output')
xxg2(1,:)



% compute  first principal component of  extracted  cycles
XX1= standard(xxo1);
OPTS.disp = 0;
[ay, by] = eigs(cov(XX1),3,'LM',OPTS);
PC1 = XX1*ay(:,1);  % the first principal component
lam1 = ay*sqrt(by);  % loadings on the PC of the variables

XX2= standard(xxo2);
[acy, bcy] = eigs(cov(XX2),3,'LM',OPTS);
PC2 = XX2*acy(:,1);  % the first principal component
lam2 = acy*sqrt(bcy);  % loadings on the PC of the variables

dates=1999.50:0.25:2017.75;
plot(dates,PC1,'r', 'linewidth',2); hold on;
plot(dates,PC2, 'k-.', 'linewidth',2);hold  off; axis  tight;
legend('PC cyclical GDP','PC cyclical Cr/GDP')
title('Real and  financial cycles')
pause
close all;

% statistics  of  PCs
xxo3=[PC1 PC2];
auto3=autocor(xxo3,4);

disp('contemporaneous correlations of  PC')
auto3(1:2,1:2)
disp('first lagged correlations of  PC')
auto3(3:4,1:2)
disp('second lagged correlations of PC')
auto3(5:6,1:2)

xxg3=[var(PC1) var(PC2)];
disp('variabilities of PC')
xxg3(1,:)

% computing  spectra
[wwm1, f1]=pwelch(PC1);
[wwm2, ~]=pwelch(PC2);

% transforming   radiant  in cycles to  position  vline
g=length(f1); p=zeros(g,1);
for  jj=1:g
    %p(jj)=fix(2*f1(g)/f1(jj));
    p(jj)=2*f1(g)/f1(jj);
end

plot(f1,wwm1,'r-', 'linewidth',2); hold on;
plot(f1, wwm2,'b--','linewidth',2); hold off; axis tight;
legend('PC1 GDP', 'PC1 credit/GDP')
title('Spectra')
vline(f1(10),'k-','p=32'); vline(f1(34),'k-','p=8');
vline(f1(5),'k-','p=64');
xlabel('Frequency')
pause
close all;



return


