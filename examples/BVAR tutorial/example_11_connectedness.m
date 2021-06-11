%% BVAR tutorial: Connectedness in the crypto currency market
% Author:   Filippo Ferroni and Fabio Canova
% Date:     4/6/2021

close all; clc;
clear all;

addpath ../../cmintools/
addpath ../../bvartools/

%% Load the data

[a,b,c] = xlsread('./crypto_all','price');
[high,~,~] = xlsread('./crypto_all','high');
[low,~,~] = xlsread('./crypto_all','low');

% define time
time = [1: size(b,1)-1]';
% time    = [datenum(b(2,1)):1:datenum(b(end,1))]';
% timestr = datestr(time);
T       = size(time,1);
timestr = b(2:end,1);


% crypto names
cryptonames = b(1,2:end);

% 1. medium cross-section
% use the 12 crypto that start in 2018/01/01 
sel1 = {'Binance Coin','Bitcoin Cash','Bitcoin','Cardano','Dash',...
'EOS','Ethereum','Ethereum Classic','IOTA','Litecoin',...
    'Monero','NEM','Neo','Qtum','Stellar','TRON','VeChain','XRP','Zcash'};
[~, indx1] = ismember(sel1,cryptonames);

% 2. larger cross-section starting in 2020/01/01
sel2 = {'Binance Coin',...
    'Bitcoin Cash',...
    'Bitcoin',...
    'Bitcoin SV',...
    'BitTorrent',...
    'Cardano',...
    'Chainlink',...
    'Cosmos',...
    'Dai',...
    'Dash',...
    'Decred',...
    'Dogecoin',...
    'EOS',...
    'Ethereum Classic',...
    'Ethereum',...
    'FTX Token',...
    'Huobi Token',...
    'IOTA',...
    'Litecoin',...
    'Maker',...
    'Monero',...
    'NEM',...
    'Neo',...
    'Polygon',...
    'Qtum',...
    'Stellar',...
    'THETA',...
    'TRON',...
    'UNUS SED LEO',...
    'USD Coin',...
    'VeChain',...
    'XRP',...
    'Zcash'};

[~, indx2] = ismember(sel2,cryptonames);
% t1 = strmatch('20-Jan-2020',timestr);
t1 = find(strcmp('1/20/2020',timestr));

% crypto log price
y1 = log(a(:,indx1));
y2 = log(a(t1:end,indx2));


%%
% remove zeros and negative prices (if any)
aa = low(1:end,indx1);
aa(aa==0) = 0.000001;
aa(aa<0 ) = 0.000001;

vol1 = 100*sqrt(365*0.361*(log(high(1:end,indx1))-log(aa)).^2);
% vol1 = 100*sqrt(365*0.361*(log(high(1:end,indx1))-log(low(1:end,indx1))).^2);
vol2 = 100*sqrt(365*0.361*(log(high(t1:end,indx2))-log(low(t1:end,indx2))).^2);

% check for nan/zeros
if any(any(isnan(y1)))
    error('there are nans in y1');
end
if any(any(isnan(y2)))
    error('there are nans in y2');
end
if any(any(isnan(vol1))) || any(any(isnan(vol2)))
    error('there are nans in vol');
end

%% plot crypto

step_plot = 300;
tmp_str= b(2:step_plot:end,1);

figure,
plot(time,log(a(:,indx2)))
set(gca,'Xtick',time(1:step_plot:end))
set(gca,'Xticklabel',tmp_str)
axis tight
title('Daily log price')
set(gcf,'position' ,[50 50 900 550])

% STR_RECAP = ['./logprice'];
% saveas(gcf,STR_RECAP,'fig');
% savefigure_pdf([STR_RECAP '.pdf']);


figure,
vol3 = real(sqrt(365*0.361*(log(high(1:end,indx2))-log(low(1:end,indx2))).^2));
plot(time, vol3)
set(gca,'Xtick',time(1:step_plot:end))
set(gca,'Xticklabel',tmp_str)
title('Daily volatility')
axis tight
ylim([-2 50])
set(gcf,'position' ,[50 50 900 550])

% STR_RECAP = ['./vol'];
% saveas(gcf,STR_RECAP,'fig');
% savefigure_pdf([STR_RECAP '.pdf']);
 

%% select the data to use for connectedness

sample = 1;% large T short N prices
% sample = 2;% short T large N prices
% sample = 3;% large T short N vol
% sample = 4;% short T large N vol

% horizon to compute the FEVD (connectedness)
options.nethor        = 10;
% options.connectedness = 2;

if sample == 1 % large T short N
    y = y1; % long times series small cross section
    cryptonames_ = cryptonames(indx1);
    tstart = 1;
    T_ = T;
    dirname = ['.\results_h' num2str(options.nethor) '\prices\largeTsmallN'];
    selection = sel1;

elseif sample ==2 % short T large N
    y = y2; % shorter times series large cross section
    cryptonames_ = cryptonames(indx2);
    tstart = t1;
    T_ = size(y2,1);
    dirname = ['.\results_h' num2str(options.nethor) '\prices\smallTlargeN'];
    %dirname = '.\results\prices\smallTlargeN';
    selection = sel2;
    
elseif sample ==3 % volatility sample 1
    y = vol1;
    cryptonames_ = cryptonames(indx1);
    tstart = 1;
    T_ = size(vol1,1);
    dirname = ['.\results_h' num2str(options.nethor) '\vols\largeTmallN'];
    %dirname = '.\results\vols\largeTsmallN';
    selection = sel1;

elseif sample ==4
%     y = vol2;
    y = log(vol2);
    cryptonames_ = cryptonames(indx2);
    tstart = t1;
    T_ = size(vol2,1);
    dirname = ['.\results_h' num2str(options.nethor) '\vols\smallTlargeN'];
%    dirname = '.\results\vols\smallTlargeN';
    selection = sel2;
end

mkdir(dirname)


%% Estimate the VAR on the full sample

% activate the Ridge estimator
options.Ridge.est         = 1;
% penalty parameter
options.Ridge.lambda      = 0.25;

% activate the Lasso estimator
options.Lasso.est         = 1;
% penalty parameter
options.Lasso.lambda      = 0.25;

% activate the ENet estimator
options.ElasticNet.est    = 1;
% penalty parameters
options.ElasticNet.lambda = 0.25;
options.ElasticNet.alpha  = 0.5;

% activate the Minnesota Prior
options.prior.name = 'Minnesota';
options.max_minn_hyper  = 1;   % start the  optimization routine
options.index_est       = 1:2;   % hyper-parameters to maximize
options.minn_prior_tau     = 10;
options.minn_prior_lambda  = 0.1;
% options.minn_prior_muu     = 0.1;

% activate the Connectedness measures
options.connectedness  = 1; % Pesaran-Shin Identification 
% options.connectedness  = 2; % Alternative Identification - default recursive
% options.signs{1}     = 'y(3,1:3,1)>0'; % BITCOIN up in period 1 to 3
% options.signs{2}     = 'y(7,1:3,1)<0'; % ETHEREUM down in period 1 to 3 

lags          = 3;
options.K     = 5000;
bvar1         = bvar_(y,lags,options);

%% print the results for the full sample

delete([dirname '\results.log'])
diary([dirname '\results.log'])
tmp0 = [mean(bvar1.Connectedness.Index) bvar1.Ridge.Connectedness.Index ...
    bvar1.Lasso.Connectedness.Index bvar1.ElasticNet.Connectedness.Index ];

tmp1 =[mean(bvar1.Connectedness.FromAlltoUnit,2) ...
    bvar1.Ridge.Connectedness.FromAllToUnit ...
    bvar1.Lasso.Connectedness.FromAllToUnit ...
    bvar1.ElasticNet.Connectedness.FromAllToUnit];

tmp2 =[mean(bvar1.Connectedness.FromUnitToAll,2) ...
    bvar1.Ridge.Connectedness.FromUnitToAll ...
    bvar1.Lasso.Connectedness.FromUnitToAll ...
    bvar1.ElasticNet.Connectedness.FromUnitToAll];

disp('============================================')
disp('Overall Connectedness Index')
disp([{'BVAR-Minnesota','Ridge','Lasso','ElasticNet'}; num2cell(tmp0)])
disp('')
disp('============================================')
disp('From Unit to All (Spillover)')
disp([{'Crypto','BVAR-Minnesota','Ridge','Lasso','ElasticNet'}; [cryptonames_' num2cell(tmp2)]])
disp('')
disp('============================================')
disp('From All to Unit (Vulnerability)')
disp([{'Crypto','BVAR-Minnesota','Ridge','Lasso','ElasticNet'}; [cryptonames_' num2cell(tmp1)]])
disp('')

diary off

%% Rolling Estimate of connectedness index

opts.Ridge.est         = 1;
opts.Ridge.lambda      = 0.5;
opts.Lasso.est         = 1;
opts.Lasso.lambda       = 0.5;
opts.ElasticNet.est    = 1;
opts.ElasticNet.lambda = 0.5;
opts.ElasticNet.alpha  = 0.5;

opts.connectedness = 1;
opts.nethor        = 10;



W       = 200;
opts.K  = 1;
dd_     = 0;
inx = nan(T_,3);
inxfromall2unit = nan(T_,length(cryptonames_),3);
inxfromunit2all = nan(T_,length(cryptonames_),3);
wb = waitbar(0, 'Rolling index');

while dd_ + W < T_
    
    dd_   = dd_+1;
    span = dd_ : W+dd_;
    bvar0         = bvar_(y(span,:),lags,opts);
    inx(dd_+W,1) = bvar0.Ridge.Connectedness.Index;
    inx(dd_+W,2) = bvar0.Lasso.Connectedness.Index;
    inx(dd_+W,3) = bvar0.ElasticNet.Connectedness.Index;
    
    inxfromall2unit(dd_+W,:,1) = bvar0.Ridge.Connectedness.FromAllToUnit;
    inxfromall2unit(dd_+W,:,2) = bvar0.Lasso.Connectedness.FromAllToUnit;
    inxfromall2unit(dd_+W,:,3) = bvar0.ElasticNet.Connectedness.FromAllToUnit;
    
    inxfromunit2all(dd_+W,:,1) = bvar0.Ridge.Connectedness.FromUnitToAll;
    inxfromunit2all(dd_+W,:,2) = bvar0.Lasso.Connectedness.FromUnitToAll;
    inxfromunit2all(dd_+W,:,3) = bvar0.ElasticNet.Connectedness.FromUnitToAll;
    
    waitbar(dd_/(T_-W), wb);
    
end

close(wb)
%%


step_plot = 200;
tmp_str = b(tstart:step_plot:end,1);

figure('Name','Rolling Connectedness Index')
plot(time(tstart:end), inx(:,1:3))
set(gca,'Xtick',time(tstart:step_plot:end))
set(gca,'Xticklabel',tmp_str)
title('Rolling Connectedness Index')
axis tight
ylim([88 99])
legend('Ridge','Lasso','ElasticNet')
STR_RECAP = [ dirname '/RolllingIndex'];
saveas(gcf,STR_RECAP,'fig');
savefigure_pdf([STR_RECAP '.pdf']);

%%

step_plot = 300;
tmp_str = b(tstart:step_plot:end,1);

for sel = 1 : length(selection)
    
    figure('Name',selection{sel})   
    
    subplot(2,2,1)
    icrypto = find(strcmp(selection{sel},cryptonames_));
    
    plot(time(tstart:end), exp(y(:, icrypto)) );
    set(gca,'Xtick',time(tstart:step_plot:end))    
    set(gca,'Xticklabel',tmp_str)
    title([selection{sel} '- USD price' ])
    axis tight
    
    subplot(2,2,2)
    plot(time(tstart:end), squeeze(inxfromunit2all(:,icrypto,1:3)))
    axis tight
    set(gca,'Xtick',time(tstart:step_plot:end))
    set(gca,'Xticklabel',tmp_str)   
    title(['From ' selection{sel}  ' to All (Spillover)'])
    
    subplot(2,2,3)
    plot(time(tstart:end), squeeze(inxfromall2unit(:,icrypto,1:3)))
    set(gca,'Xtick',time(tstart:step_plot:end))
    set(gca,'Xticklabel',tmp_str)
    title(['From All to ' selection{sel} ' (Vulnerability)'])
    axis tight
    legend('Ridge','Lasso','ElasticNet')   
    
    set(gcf,'position' ,[50 50 900 450])

    STR_RECAP = [ dirname '/RolllingIndex_' selection{sel}];
    saveas(gcf,STR_RECAP,'fig');
    savefigure_pdf([STR_RECAP '.pdf']);

    
end


