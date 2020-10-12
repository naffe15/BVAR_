clear all
close all
clc

addpath ../../cmintools/
addpath ../../v4.1


% dating_exa.m, 16-02-2020, version 1.3, fabio canova
% this program illustrates the use of BB program to  date  turning  points
% and  compute business cycle statistics

%%%%%%%%%%%%%%%%%%%%%%% exercises  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  exercise 1: Vary the censoring  rules.
%  exercise 2: Use  a  shorter  sample. Do  the  turning  points  coincide?
%  exercise 3: Change  the  thresh  parameter.
%  exercise 4: Change  the  data  set.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=======================================================================
% LOAD DATA
%=======================================================================
lload=0; % =0 euro data, =1 US data

if lload==0
    % Euro area AWM DATABASE: Quarterly
    [a,b,~] = xlsread('awm19up18.csv');
    % names of variables
    varnames = b(1,2:end);
    
    % time convention: Q1 = .00 and Q4 =0.75
    time = 1970 : .25 : 2017.75;
    time_start   = find(time==1970.50);
    time_start1  = find(time==1999.50);
    time_end     = find(time==2017.75);
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
    
    %  data for  credit  to  GDP starts  only  at  time_start1 (1999:25)
    ddata(:,11) = yy(time_start+1:time_end,11);
    endd=length(ddata);
    
    % pick log output, log consumption, log  investment, log labor
    % productivity, log employment, interest rate, inflation rate
    ddd=ddata(:,[1 2 3 9 10 6 5]);
    T=length(ddd);
    time_data = time(time_start+1:time_end);
    
    
else
    
    % US DATABASE: Quarterly
    
    % real GDP, urate, real consumption, gdp defl, real  investment,
    % capacity utilization, call rate, 10y goverment  bond  rate
    % names: RGDP, URATE, C, GDPdef, inv, capU, callrate, 10ygbond
    [c,d,~] = xlsread('USdata.xlsx');
    % names of variables
    varnames = d(1,2:end);
    
    % time convention: Q1 = .00 and Q4 =0.75
    time = 1970 : .25 : 2019.50;
    time_start   = find(time==1969.75);
    time_end     = find(time==2017.75);
    time_break   = find(time==2007.75);
    
    tt=length(c);
    
    % data transformations:
    
    % 1. difference log of output, consumption,investment
    yy1(1: tt,[1 2 3]) = log(c(1 : tt,[1 3 5]));
    % 2. difference the log GDP deflator index
    yy1(2: tt,6) = diff(log(c(1: tt,4)))*400;
    % 3. leave  the Urate capU unchanged
    yy1(2: tt,[4 5])= c(2:tt,[2 6]);
    % 4. leave call rate, 10 year unchanged
    yy1(2: tt,7:8)= c(2:tt,7:8);
    % 5. compute detrended C/Y, I/Y ratio
    yy1(1: tt,9) = c(1 : tt,3)./c(1: tt,1);
    yy1(1: tt,10) = c(1 : tt,5)./c(1: tt,1);
    % 6. Compute term  spread
    yy1(1: tt,11) = c(1 : tt,8)-c(1: tt,7);
    
    ddd=yy1;
    T=length(ddd);
end


%% Parameters

%frequency
freq   = 'q';        % 'q' for quarterly, 'm' for monthly %
tstart = [1970 3];
tend   = [2017 3];

% cycle parameters
options.turnphase   = 2;
options.phase       = 2;          % censoring rules %
options.cycle       = 5;          % lenght of cycle
options.thresh      = 10.4;       % bypasses phase and cycle restriction if peak to trough is > than thresh

options.nrep     = 1;        % 1 if analyze real data
options.complete = 1;        % if= 1- use complete cycles,if =0 -use incomplete cycles (excess still on complete cycle)


dura=zeros(size(ddd,2),2);  ampl=zeros(size(ddd,2),2);
cumm=zeros(size(ddd,2),2);  excc=zeros(size(ddd,2),2);
durcv=zeros(size(ddd,2),2); ampcv=zeros(size(ddd,2),2);
exccv=zeros(size(ddd,2),2); nott=zeros(size(ddd,2),1);
turn    = time(time_start+1:time_end)';
zz      = time(time_start+1:time_end)';
tid     = linspace(1970.75,2017.75,T)';

for qq = 1 : size(ddd,2)
    
    x = ddd(:,qq);
    
    % dating_mbbq
    disp('series')
    disp(qq)

    [dt_(qq)] = date_(x, time_data, freq, tstart, tend, options); 
        
    zz      = [zz dt_.st];
    turn    = [turn dt_.trinary];
    
    dura(qq,:)  = dt_.dura;
    ampl(qq,:)  = dt_.ampl;
    cumm(qq,:)  = dt_.cumm;
    excc(qq,:)  = dt_.excc;
    durcv(qq,:) = dt_.durcv;
    ampcv(qq,:) = dt_.amplcv;
    exccv(qq,:) = dt_.exccv;
    nott(qq,:)  = dt_.notentp;
    %plot(tid,trinary,'linewidth',1)
    %title('Turning  points')
    
end

disp('statistics on average cycle')

disp('duration contractions/duration expansions')
disp(dura)

disp('amplitudes contractions/amplitude expansions')
disp(ampl)

disp('cumulative contractions/cumulative expansions')
disp(cumm)

disp('excess movements percent of triangle area')
disp('contractions/expansions')
disp(excc)

disp('cv of durations contractions/expansions')
disp(durcv)

disp('cv of amplitudes contractions/expansions')
disp(ampcv)

disp('cv of excess movements contractions/expansions')
disp(exccv)

disp('no of its skipped since no peaks+troughs<=2')
disp(nott)

%disp('states indicators:contraction=0, expansion=1')
%format  short  g
%round(zz, 2)

disp('concordance index BC phases')
ma=corr(zz,'type','Spearman');
disp(ma(2,3:size(ma,2)))

disp('concordance index turning  points')
qa=corr(turn,'type','Spearman');
disp(qa(2,3:size(qa,2)))


[aa,bb]=size(turn);
tsum=zeros(aa,1);
for kk=1:aa
    ssum=0;
    for  qq=2:bb
        ssum=ssum+turn(kk,qq);
    end
    tsum(kk)=ssum;
end

plot(tid,tsum,'linewidth',2)
title('Distribution of Turning  points')

disp('Distribution of  peaks')
dp=[zz(find(tsum>=1),1) tsum(find(tsum>=1),1)];
disp(dp)
disp('Distribution of  throughs')
dt=[zz(find(tsum<=-1),1) tsum(find(tsum<=-1),1)];
disp(dt)

% dates
% peaks    mean 1974.00; 1980.00; 1992.00; 2001.00;  2008:00;  2011.25;
%          mode 1974.50; 1980.00; 1992.00; 2001.25;  2008.00;  2011.50;
% throughs mean 1975.25; 1984.25; 1993.75; 2002.75;  2009:50;  2013.00;
%          mode 1975.00; 1984.25; 1993.75; NaN;      2010.00;  2013.00;


% costructing a  recession indicator
recind=zeros(time_end-time_start,1);

% recession dates
rec1b = find(time==1974.00);
rec1e = find(time==1975.25);

rec2b = find(time==1980.00);
rec2e = find(time==1984.25);

rec3b = find(time==1992.00);
rec3e = find(time==1993.75);

rec4b = find(time==2001.00);
rec4e = find(time==2002.75);

rec5b = find(time==2008.00);
rec5e = find(time==2009.50);

rec6b = find(time==2011.25);
rec6e = find(time==2013.00);


for i=rec1b:rec1e
    recind(i,1)=1;
end

for i=rec2b:rec2e
    recind(i,1)=1;
end
for i=rec3b:rec3e
    recind(i,1)=1;
end

for i=rec4b:rec4e
    recind(i,1)=1;
end

for i=rec5b:rec5e
    recind(i,1)=1;
end

for i=rec6b:rec6e
    recind(i,1)=1;
end

if  lload==0
    save Eurorec recind
else
    save Usarec recind
end
return









