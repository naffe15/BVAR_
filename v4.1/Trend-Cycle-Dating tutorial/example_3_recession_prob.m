clear all
close all
clc

addpath ../../cmintools/
addpath ../../v4.1

% example_3_recession_probability.m, Fabio Canova, 15/03/2020
% predicting  recession  and  the start of  a recession 
% in the  euro  area. use recession  indicator  constructed  with
% example_2_dating.m and  a few  regressors ( labor  productivity, commodity
% prices, long term  rates  or spread).


% just  checking  that  the  program  works
%y=rand(100,1); x=zeros(100,1);
%for  i=1:100
%    if y(i)<0.5
%        x(i)=0;
%    else
%        x(i)=1;
%    end
%end
%z =randn(100,3);
%[estimator,cov_Hessian,ME1,ME2,ME_std] = probit(x,z);


% Euro area AWM DATABASE: Quarterly
    [a,b,~] = xlsread('awm19up18.csv');
    % names of variables
    varnames = b(1,2:end);

    % time convention: Q1 = .00 and Q4 =0.75
    time = 1970 : .25 : 2017.75;
    time_start = find(time==1970.50); 
    time_start1 = find(time==1999.50); 
    time_end   = find(time==2017.75);
    time_break   = find(time==2007.75);
  
    
    % The  CREDIT DATA: Quartely
    [e,f,~]=xlsread('ECB_Credit_gaps.xlsx');
    varnames2= f(1,2:end);

    % real GDP, real consumption, real  investment, GDP defl, HICP, 
    % short and long term interest rate, commodity  prices, 
    % labor  productivity, total  employment, credit to NFC to  GDP 
    yy = [a(:,1) a(:,2) a(:,4) a(:,7) a(:,19) a(:,33) a(:,34) ...
        a(:,35) a(:,42) a(:,29) e(time_start-2:time_end,25)]; 
    % names: YER PCR ITR YED HICP STN LTN COMPR LPROD LNN

    % data transformations:
    % 1. the log of output, consumption,investment
    ddata(:,1:3) = log(yy(time_start+1 : time_end,1:3));
    % 2. the log/log difference of GDP deflator index and CPI
    %ddata(:,4:5) = log(yy(time_start+1  : time_end,4:5));
    ddata(:,4:5) = diff(log(yy(time_start : time_end,4:5)))*400;
    % 3. the level of short and long term interest rate
    ddata(:,6:7) = yy(time_start +1 : time_end,6:7);
    % 4. the log/log difference of commodity prices
    ddata(:,8) = log(yy(time_start+1  : time_end,8));
    %ddata(:,8) = diff(log(yy(time_start : time_end,8)))*400;
    % 5. the taking the log of labor  productivity and  employment
    ddata(:,9:10) = log(yy(time_start+1 : time_end,9:10));
    %ddata(:,9) = diff(log(yy(time_start+1 : time_end,9)));
    %ddata(:,10) = log(yy(time_start+1 : time_end,10));

    % 6. spread
    ddata(:,11)= ddata(:,7)-ddata(:,6);
    %  data for  credit  to  GDP starts  only  at  time_start1 (1999:25)
    ddata(:,12) = yy(time_start+1:time_end,11);
    endd=length(ddata); 

 % cc=ones(time_end-time_start,1);  % constant
 % pick log labor productivity, log commodity prices, long term rate
 % z=ddata(:,[9 8 7]);
 % pick log labor productivity, log commodity prices, spread
 z=ddata(:, [9 8 11]);
    
load  Eurorec
% recind has  the  recession indicator  created  with  dating_exa.m
x=recind;

% use contemporanous  values
% [estimator,cov_Hessian,ME1,ME2,ME_std] = probit(x,z);

% use  one  period lagged  values  
zz=zeros(size(z,1), size(z,2));
zz(2:length(z),:)=z(1:length(z)-1,:);
timeplot = time(time_start:time_end-1);

% Estimate Probit model and the marginal effects
[estimator,cov_Hessian,ME1,ME2,ME1_std] = probit(x,zz);

% prediction
pz=[ones(length(zz),1) zz];
px=normpdf(pz*estimator);
half=0.5*ones(length(z),1);

plot(timeplot,px,'r-.','Linewidth',2); hold on; 
plot(timeplot,half,'b:','Linewidth',2); hold  on; 
plot(timeplot,x,'k-','Linewidth',2); hold  off; axis  tight;
legend('predition','halfline', 'actual')
pause
close  all;

% do  estimation recursively: predicting  the  beginning  of  a  recession
% which is  the  same  as  predicting  a  peak
% recession dates
%rec1b = find(time==1974.00);
%rec1e = find(time==1975.25);
%[estimator1,cov_Hessian1,ME11,ME21,ME11_std] = probit(x(1:rec1b),zz(1:rec1b,:));

rec2b = find(time==1980.00);
rec2e = find(time==1984.25);
disp('Predicting 1980 recession')
[estimator2,cov_Hessian2,ME12,ME22,ME12_std] = probit(x(1:rec2b),zz(1:rec2b,:));
pz2=[ones(rec2b,1) zz(1:rec2b,:)];
px2=normpdf(pz2*estimator2);
half2=0.5*ones(rec2b,1);

plot(timeplot(1:rec2b),x(1:rec2b),'k-','Linewidth',2); hold  on; 
plot(timeplot(1:rec2b),half2(1:rec2b),'b:','Linewidth',2); hold  on; 
plot(timeplot(1:rec2b),px2(1:rec2b),'r-.','Linewidth',2); hold  off;axis tight;
legend('actual', 'halfline', 'predition')
pause
close all;

rec3b = find(time==1992.00);
rec3e = find(time==1993.75);
disp('Predicting  1992 recession')
[estimator3,cov_Hessian3,ME13,ME23,ME13_std] = probit(x(1:rec3b),zz(1:rec3b,:));
pz3=[ones(rec3b,1) zz(1:rec3b,:)];
px3=normpdf(pz3*estimator3);
half3=0.5*ones(rec3b,1);

plot(timeplot(1:rec3b),x(1:rec3b),'k-','Linewidth',2); hold  on; 
plot(timeplot(1:rec3b),half3(1:rec3b),'b:','Linewidth',2); hold  on; 
plot(timeplot(1:rec3b),px3(1:rec3b),'r-.','Linewidth',2); hold  off; axis  tight;
legend('actual', 'halfline','predition')
pause
close all;

rec4b = find(time==2001.00);
rec4e = find(time==2002.75);
disp('Predicting  2001 recession')
[estimator4,cov_Hessian4,ME14,ME24,ME14_std] = probit(x(1:rec4b),zz(1:rec4b,:));
pz4=[ones(rec4b,1) zz(1:rec4b,:)];
px4=normpdf(pz4*estimator4);
half4=0.5*ones(rec4b,1);

plot(timeplot(1:rec4b),x(1:rec4b),'k-','Linewidth',2); hold  on; 
plot(timeplot(1:rec4b),half4(1:rec4b),'b:','Linewidth',2); hold  on; 
plot(timeplot(1:rec4b),px4(1:rec4b),'r-.','Linewidth',2); hold  off; axis  tight;
legend('actual', 'halfline','predition')
pause
close all;


rec5b = find(time==2008.00);
rec5e = find(time==2009.50);
disp('Predicting 2008 recession')
[estimator5,cov_Hessian5,ME15,ME25,ME15_std] = probit(x(1:rec5b),zz(1:rec5b,:));
pz5=[ones(rec5b,1) zz(1:rec5b,:)];
px5=normpdf(pz5*estimator5);
half5=0.5*ones(rec5b,1);

plot(timeplot(1:rec5b),x(1:rec5b),'k-','Linewidth',2); hold  on; 
plot(timeplot(1:rec5b),half5(1:rec5b),'b:','Linewidth',2); hold  on; 
plot(timeplot(1:rec5b),px5(1:rec5b),'r-.','Linewidth',2); hold  off; axis tight;
legend('actual', 'halfline', 'predition')
pause
close all;


rec6b = find(time==2011.25);
rec6e = find(time==2013.00);
disp('Predicting  2011 recession')
[estimator6,cov_Hessian6,ME16,ME26,ME16_std] = probit(x(1:rec6b),zz(1:rec6b,:));
pz6=[ones(rec6b,1) zz(1:rec6b,:)];
px6=normpdf(pz6*estimator6);
half6=0.5*ones(rec6b,1);

plot(timeplot(1:rec6b),x(1:rec6b),'k-','Linewidth',2); hold  on; 
plot(timeplot(1:rec6b),half6(1:rec6b),'b:','Linewidth',2); hold  on; 
plot(timeplot(1:rec6b),px6(1:rec6b),'r-.','Linewidth',2); hold  off; axis  tight;
legend('actual', 'halfline', 'predition')
 
 
return
    
    
 