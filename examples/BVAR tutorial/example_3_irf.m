%% BVAR tutorial: Impulse  responses, variance  decomp, historical decomp
%                 plus  time varying  IRFs
% Authors:   Filippo Ferroni and  Fabio Canova
% Date:     27/02/2020, revised  14/12/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Computation of IRFs using various identification schemes
% 2) Calculation of  variance and  historical decompositions
% 3) Rolling window estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clc; clear;
addpath ../../cmintools/
addpath ../../bvartools/

% load the data
load DataGK
y = [logip logcpi gs1 ebp];

%% 1) Cholesky identification

% run the BVAR
lags        = 12;
options.hor = 24;
bvar1       = bvar_(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho              = 3;   
% Change the order of the variables for the plot 
% 1. gs1; 2. logcpi; 3. logip; 4. ebp
indx_var              = [3, 2, 1, 4];
% IRF to PLOT
irfs_to_plot           = bvar1.ir_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% variables names for the plots
options.varnames      = {'1 year rate','CPI','IP','EBP'};  
% names of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% names of the figure to save
options.saveas_strng  = 'cholesky';
% name of the shock
options.shocksnames   = {'MP'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;   
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)
pause;

%% 2) Sign restrictions

% specify the restrictions
options.signs{1} = 'y(3,1:3,1)>0'; % 1Y rate up in period 1 to 3
options.signs{2} = 'y(2,1:3,1)<0'; % CPI down in period 1 to 3 
options.K        = 1000;
% run the BVAR
bvar2             = bvar_(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho         = 1;  
% IRF  (notice that we change the order of the variables for the plot)
irfs_to_plot     = bvar2.irsign_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% names of the figure to save
options.saveas_strng  = 'signs';
% name of the shock
options.shocksnames   = {'MPtightening'}; % 
plot_irfs_(irfs_to_plot,options)
pause;


%% 3) Signs and Narrative restrictions

% specify the restrictions
options.signs{1}     = 'y(3,1:3,1)>0'; % 1Y rate up in period 1 to 3
options.signs{2}     = 'y(2,1:3,1)<0'; % CPI down in period 1 to 3 
% Large Tighening episodes of Volker are 
% m11-1980	1.871416916
% m05-1981	1.514907827
% see Narrative Romer and Romer
% Sample starts in 1979m7 and there are 12 lags. The first reduced form 
% innovation is 1980m7. Therefore, 
% 1980m09-11 = # 3:5
% 1981m05 = # 11
options.narrative{1} = 'v([3:5],1)>0'; 
options.narrative{2} = 'v([11],1)>0'; 

% index = find(T==2007+8/12)-lags; % september 2007
% options.narrative{1} = ['v([ ' num2str(index)  ' ],1)<0']; 

% run the BVAR
options.K        = 1000;
bvar3             = bvar_(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho         = 1;  
% IRF (notice that we change the order of the variables for the plot)
irfs_to_plot     = bvar3.irnarrsign_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% names of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% names of the figure to save
options.saveas_strng  = 'signsnarrative';
% name of the shock
options.shocksnames   = {'MPtightening'}; % 
plot_irfs_(irfs_to_plot,options)
pause;

%% 4) Zeros and Signs

% Housekeeping: remove previous identification settings
options = rmfield(options,'signs');
options = rmfield(options,'narrative');

% specify the restrictions
% 1) ad = aggregate demand disturbance [sign restrictions]
options.zeros_signs{1}     = 'y(1,1)=1;';
options.zeros_signs{end+1} = 'y(2,1)=1;'; 
options.zeros_signs{end+1} = 'y(3,1)=1;';
% 2) as = aggregate supply shock [sign restrictions]
options.zeros_signs{end+1} = 'y(1,2)=1;';
options.zeros_signs{end+1} = 'y(2,2)=-1;';
% 3) mp = monetary policy shock, no cont response of prices and quantities 
% [zero restrictions]
options.zeros_signs{end+1} = 'ys(1,3)= 0;';
options.zeros_signs{end+1} = 'ys(2,3)=0;';
% 3) mp = rate and bond premium go up [sign restrictions]
options.zeros_signs{end+1} = 'y(3,3)=1;';
options.zeros_signs{end+1} = 'y(4,3)=1;';
% run the BVAR
bvar4             = bvar_(y,lags,options);

% IRFs of Interest
indx_sho     = 1:3;
irfs_to_plot = bvar4.irzerosign_draws(indx_var,:,indx_sho,:);

% Customize the IRF plots
options.saveas_strng    = 'zerossigns';
options.shocksnames     = {'ADshock','ASshock','MPshock'}; % 
% options.add_irfs        = bvar4.irzerosign_ols(indx_var,:,indx_sho,:);
plot_all_irfs_(irfs_to_plot,options)
pause

%% 5) Long Run Restrictions: Technology shock

% Housekeeping: remove previous identification settings
options = rmfield(options,'zeros_signs');
options = rmfield(options,'saveas_strng');
options = rmfield(options,'shocksnames');
options.K        = 1000;

% define the dataset for the identification of LR shock
y = [diff(logip) logcpi(2:end) gs1(2:end) ebp(2:end)];

% run the BVAR
options.long_run_irf   = 1; 
% options.zeros_signs{1} = 'yl(1,1)=0;';
bvar5                  = bvar_(y,lags,options);

% IRF of Interest
% shock index
indx_sho     =  1;
% Define the order of the variables for the plot 
% 1. D(logip); 2. logcpi; 3.  gs1; 4. ebp;
indx_var     = 1:4;
irfs_to_plot          = bvar5.irlr_draws(indx_var,:,indx_sho,:);
cirfs_to_plot         = cumsum(bvar5.irlr_draws(indx_var,:,indx_sho,:),2);
irfs_to_plot(1,:,:,:) = cirfs_to_plot(1,:,:,:);

% Customize the IRF plot
% variables names for the plots
options.varnames      = {'IP','CPI','1 year rate','EBP'};  
% names of the figure to save
options.saveas_strng    = 'LR';
% name of the shock
options.shocksnames     = {'Technology'}; % 
plot_irfs_(irfs_to_plot,options)
pause;


%% 6) Proxy or IV 

% define the dataset for the identification of MP shock with IV
y = [gs1 logip logcpi ebp];

% load the instruments
datafile = 'factor_data.xlsx';
sheet    = 'factor_data';
if isMATLABReleaseOlderThan("R2024a")
    [numi,txti,rawi] = xlsread(datafile,sheet);
    % use the same instrument as GK
    options.proxy  = numi(:,4);
else
    Tbl = readtable(datafile,Sheet=sheet);
    % use the same instrument as GK
    options.proxy  = Tbl.ff4_tc;
end

% run the BVAR
bvar6        = bvar_(y,lags,options);

% IRF of Interest
% shock index
indx_sho     =  1;
% Keep the same order of variables as in the estimation
% 1. gs1; 2. logip; 2. logcpi; 4. ebp;
indx_var     = 1:4;
irfs_to_plot = bvar6.irproxy_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% variables names for the plots
options.varnames      = {'1 year rate','IP','CPI','EBP'};  
% names of the figure to save
options.saveas_strng    = 'IV';
% name of the shock
options.shocksnames     = {'MP'}; 
plot_irfs_(irfs_to_plot,options)

% %%
% 
% % consider the mean of the posterior distribution 
% Phi   = mean(bvar1.Phi_draws,3);
% Sigma = mean(bvar1.Sigma_draws,3);
% i =1; j=1;
% 
% crit = nan(100,1);
% Q    = nan(bvar1.N,bvar1.N,100);
% for k = 1 : 100
%     Q(:,:,k)  = generateQ(bvar1.N);               % generate an orthonormal matrix
%     FEVD     = fevd(h,Phi,Sigma,Q(:,:,k));  % compute the FEVD
%     % Contribution of shock j to variable i forecast error volatility at horizon h
%     crit(k,1) = FEVD(i, h, j);
% end
% % Pick the argmax
% [~,index] = max(crit);
% Qbar      = Q(:,:,index);



%% Extra part 1): Forecast Error variance decomposition with Cholesky (bvar1)
% compute the contribution of MP to the H-step ahead forecast error  of
% observables

% consider the mean of the posterior distribution 
Phi   = mean(bvar1.Phi_draws,3);
Sigma = mean(bvar1.Sigma_draws,3);
% index of the shocks of interest (shock to gs1)
indx_sho              = 3;   
% 2 year ahead FEVD
hh      = 24;
FEVD    = fevd(hh,Phi,Sigma);

disp('%=====================================================%')
disp('% Forecast Error Variance Decomposition               %')
disp('% Percentage of volatility explained by MP shock      %')
disp('    logip     logcpi   gs1        ebp ')
disp(FEVD(:,indx_sho)')
disp('%                                                     %')
disp('%=====================================================%')
pause;


%% Extra part 2): Historical Decomposition 
% VAR with zero-sign restrictions
clear
load DataGK
y = [logip logcpi gs1 ebp];
lags        = 12;


% specify the restrictions
% 1) ad = aggregate demand disturbance [sign restrictions]
options.zeros_signs{1}     = 'y(1,1)=1;';
options.zeros_signs{end+1} = 'y(2,1)=1;'; 
options.zeros_signs{end+1} = 'y(3,1)=1;';
% 2) as = aggregate supply shock [sign restrictions]
options.zeros_signs{end+1} = 'y(1,2)=1;';
options.zeros_signs{end+1} = 'y(2,2)=-1;';
% 3) mp = monetary policy shock, no cont response of prices and quantities 
% [zero restrictions]
options.zeros_signs{end+1} = 'ys(1,3)= 0;';
options.zeros_signs{end+1} = 'ys(2,3)=0;';
% 3) mp = rate and bond premium go up [sign restrictions]
options.zeros_signs{end+1} = 'y(3,3)=1;';
options.zeros_signs{end+1} = 'y(4,3)=1;';
% run the BVAR
option.k=1000;
bvar4             = bvar_(y,lags,options);



% uses the first accepted rotation for  illustration
opts_.Omega         =  bvar4.Omegaz(:,:,1); 
[yDecomp,ierror]  = histdecomp(bvar4,opts_); 
% yDecomp = historical decomposition
% time, variable, shocks (& exogenous variables if any) and initial
% condition (last)
% ierror = structural innovation



% Declare the names of the variables in the order they appear in the VAR
bvar4.varnames      = {'IP','CPI','Interest Rate','EBP'};
% select the variables for the plot of the historical decomposition
optnsplt.plotvar_    = {'Interest Rate','EBP'};
% select the shocks combination to report
optnsplt.snames_ = { {'Shck1','Shck2'};...    Combine Supply and Demand
    {'Shck3'};...              MP    
    {'Shck4'} ...              Other shock not identified in the VAR
    };
% declare the name of the shocks for the legend
optnsplt.stag_       = {'Supply+Demand shocks';
            'MP shocks';
            'Other Shocks';
            'Deterministic component'};
% name of the file to save
optnsplt.save_strng    = 'y0';
% define the time for the plot
optnsplt.time          = T(1+lags:end);
% define the directory where the plot is saved 
optnsplt.saveas_dir    = './sdcmp_plt';
% limit the plot to a specific time window 
optnsplt.Tlim          = [2006 2012];
plot_sdcmp_(yDecomp,bvar4,optnsplt)
pause;


% loop over accepted rotations
% WARNING:  IT  TAKES  A  SEVERAL  MINUTES TO RUN
%{
    for  ii=1:100
        opts_.Omega         =  bvar4.Omegaz(:,:,ii);      
        [yDecompl,ierror]  = histdecomp(bvar4,opts_);         
        % yDecomp = historical decomposition
        % time, variable, shocks (& exogenous variables if any) and initial
        % condition (last)
        % ierror = structural innovation        
        yDecompb(:,:,:,:,ii)=yDecompl(:,:,:,:);
    end
    yDecomp=median(yDecompb,4);

% Declare the names of the variables in the order they appear in the VAR
bvar41.varnames      = {'IP','CPI','Interest Rate','EBP'};
% select the variables for the plot of the historical decomposition
optnsplt.plotvar_    = {'Interest Rate','EBP'};
% select the shocks combination to report
optnsplt.snames_ = { {'Shck1','Shck2'};...    Combine Supply and Demand
    {'Shck3'};...              MP    
    {'Shck4'} ...              Other shock not identified in the VAR
    };
% declare the name of the shocks for the legend
optnsplt.stag_       = {'Supply+Demand shocks';
            'MP shocks';
            'Other Shocks';
            'Deterministic component'};
% name of the file to save
optnsplt.save_strng    = 'y0';
% define the time for the plot
optnsplt.time          = T(1+lags:end);
% define the directory where the plot is saved 
optnsplt.saveas_dir    = './sdcmp_plt';
% limit the plot to a specific time window 
optnsplt.Tlim          = [2006 2012];
plot_sdcmp_(yDecomp,bvar41,optnsplt)
pause;
%}


%% Extra part 3): Minnesota Priors IRF 

clear;

load DataGK
y = [logip logcpi gs1 ebp];
lags                = 12;
options.hor         = 48;
bvar1               = bvar_(y,lags,options);
options.presample      = 12;
options.prior.name     = 'Minnesota';
options.minn_prior_tau = 0.5;
bvar6                  = bvar_(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho              = 3;   
% Change the order of the variables for the plot 
% 1. gs1; 2. logcpi; 3. logip; 4. ebp
indx_var              = [3, 2, 1, 4];
% IRF to PLOT
irfs_to_plot           = bvar6.ir_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% variables names for the plots
options.varnames      = {'1 year rate','CPI','IP','EBP'};  
% names of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% names of the figure to save
options.saveas_strng  = 'BayesianCholesky';
% name of the shock
options.shocksnames   = {'MP'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;
% add the Choleskyi IRF with flat prior
options.add_irfs = squeeze(median(bvar1.ir_draws(indx_var,:,indx_sho,:),4));
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)
pause;



%% EA MP: QE shocks - identification via heteroskedasticity
clear; close all; clc

% EA MP events database + Macro data
% Altavilla et al. / Journal of Monetary Economics 108 (2019) 162–179
load DataEAMP.mat
y                 = [DE10Y-DE2Y IPI HICP CORE LTR10Y];
options.varnames  = {'DE10Y-DE2Y','IPI','HICP','CORE','LRT10Y'}; 

% Identification via heteroskedasticity - pre crisis period less volatile
figure,
shade(2008+0/12,2019+11/12,[0.85 0.85 0.85],max(DE10Y-DE2Y),min(DE10Y-DE2Y))
hold on
plot(T,DE10Y-DE2Y,'k','Linewidth',2);
ylabel('DE10Y-DE2Y')
title('Changes in the long-end of the yield curve on MP events')
axis tight

% Options for Volatility-Changes identification
opts.heterosked_regimes = zeros(length(T),1);
% high volatility regime starts in 2008
opts.heterosked_regimes(find(T==2008) : end, 1) = 1;
opts.hor = 48;
lags     = 6;
bvar7    = bvar_(y,lags,opts);

% IRF to plot:
indx_var        = 1:size(y,2);
indx_sho        = 1;
irfs_to_plot    = bvar7.irheterosked_draws(indx_var,:,indx_sho,:);
% normalize to a 25 bpt increase in the 10-2Y spread
irfs_to_plot    = 0.25*irfs_to_plot/median(squeeze(irfs_to_plot(1,1,1,:)));

% names of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% names of the figure to save
options.saveas_strng  = 'Heterosk';
% name of the shock
options.shocksnames   = {'MP'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;
options.nplots        = [2 3];
% add the Cholesi IRF with flat prior
%XXX = squeeze(median(bvar7.ir_draws(indx_var,:,indx_sho,:),4));
%options.add_irfs =  0.25*XXX/XXX(1,1,1);
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)
pause;

%% Rolling window estimation 
%%
% housekeeping
clear
% load Qaurterly data
load DataQ
% collect the data: GDP deflator YoY inflation, Unemployment, 3m Tbill 
yQ          = [GDPDEF_PC1 UNRATE TB3MS];
lags        = 2;
Wsize       = 120;      % lenght of the rolling window
shift       = 4;       % time shift between adiacent windows  
w           = 0;       % index four caounting the windows 
indx_sho    = 3;       % shocks of interest
options.K   = 1;       % no draws needed we use the OLS IRF (BVAR.ir_
options.hor = 24;      % IRF horizon
rollIRF     = ...      % initialize rolling IRF (One IRF per window)
    nan(size(yQ,2),options.hor,1);
timespan       = nan(1,Wsize);

while w*shift + Wsize < size(yQ,1) 
    w = w + 1;
    timespan(w,:) = shift*(w-1) + 1 :  Wsize + shift * (w-1);
    rollbvar   = bvar_(yQ(timespan(w,:),:),lags,options);
    % normalize for the size of the shock
    norm           = rollbvar.ir_ols(indx_sho,1,indx_sho);
    % Collect the MP OLS response in the window
    rollIRF(:,:,w) =  squeeze(rollbvar.ir_ols(:,:,indx_sho))/norm;
end

figure('name','UNR')
tt = T(timespan(:,end)');
surf(tt,1:options.hor,squeeze(rollIRF(2,:,:)))
axis tight
pause;



%%% CREDIBLE SET a-la Giacomini-Kitagawa
%%% IT TAKES A  FEW  MINUTES  TO  RUN

%% Sign restrictions with Robust Credible sets 
clear; close all; clc;

% load the Gertler and Karadi data
load DataGK
y = [logip logcpi gs1 ebp];
lags = 12;

% sign restrictions
options.signs{1} = 'y(3,1:6,1)>0'; % 1Y rate up for 6 months
options.signs{2} = 'y(2,1:6,1)<0'; % CPI down for 6 months
options.K   = 1000;
options.hor = 24;
% This options activate the Kitagawa-Giacomini robust confidence set
options.robust_credible_regions.KG = 1;
% Number of draws to approximate the identified set for each posterior draw
options.robust_credible_regions.L = 100; %[default = 1000]
% The credibility level is set to .68
% options.robust_credible_regions.aalpha < 1; [default = 0.68]
% number of points on discrete grid to compute credible regions
% options.robust_credible_regions.gridLength = integer; [default = 1000]
bvar2r             = bvar_(y,lags,options);

% Define the IRF of Interest
% Change the order of the variables for the plot 
% 1. gs1; 2. logcpi; 3. logip; 4. ebp
indx_var              = [3, 2, 1, 4];
% index of the shocks of interest (shock to gs1)
indx_sho         = 1;  
% IRF  (notice that we change the order of the variables for the plot)
% with standard confidence set (a single prior)
irfs_to_plot     = bvar2r.irsign_draws(indx_var,:,indx_sho,:);
% 68% robust credible set 
options.add_irfs(:,:,1,1)  = bvar2r.irsign_robust_credible_bands.l(indx_var,:,indx_sho);
options.add_irfs(:,:,1,2)  = bvar2r.irsign_robust_credible_bands.u(indx_var,:,indx_sho);


% Customize the IRF plot
options.varnames      = {'1 year rate','CPI','IP','EBP'};  
% names of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% names of the figure to save
options.saveas_strng  = 'signs-rob';
% name of the shock
options.shocksnames   = {'MPtightening'}; % 
plot_irfs_(irfs_to_plot,options)
pause;

%% Zero-Sign restrictions with Robust Credible sets

% Housekeeping: remove previous identification settings
options = rmfield(options,'signs');
options = rmfield(options,'add_irfs');

% specify the restrictions
% 1) ad = aggregate demand disturbance [sign restrictions]
options.zeros_signs{1}     = 'y(1,1)=1;';
options.zeros_signs{end+1} = 'y(2,1)=1;'; 
options.zeros_signs{end+1} = 'y(3,1)=1;';
% 2) as = aggregate supply shock [sign restrictions]
options.zeros_signs{end+1} = 'y(1,2)=1;';
options.zeros_signs{end+1} = 'y(2,2)=-1;';
% 3) mp = monetary policy shock, no cont response of prices and quantities 
% [zero restrictions]
options.zeros_signs{end+1} = 'ys(1,3)= 0;';
options.zeros_signs{end+1} = 'ys(2,3)=0;';
% 3) mp = rate and bond premium go up [sign restrictions]
options.zeros_signs{end+1} = 'y(3,3)=1;';
options.zeros_signs{end+1} = 'y(4,3)=1;';

% This options activate the Kitagawa-Giacomini robust confidence set
options.robust_credible_regions.KG = 1;
% Number of draws to approximate the identified set for each posterior draw
options.robust_credible_regions.L = 100; %[default = 1000]
% The credibility level is set to .68
% options.robust_credible_regions.aalpha < 1; [default = 0.68]
% number of points on discrete grid to compute credible regions
% options.robust_credible_regions.gridLength = integer; [default = 1000]
% run the BVAR
bvar4r             = bvar_(y,lags,options);

% IRFs of Interest
indx_sho     = 1:3;
irfs_to_plot = bvar4r.irzerosign_draws(indx_var,:,indx_sho,:);
% 68% robust credible set 
options.add_irfs(:,:,:,1)  = bvar4r.irzerosign_robust_credible_bands.l(indx_var,:,indx_sho);
options.add_irfs(:,:,:,2)  = bvar4r.irzerosign_robust_credible_bands.u(indx_var,:,indx_sho);

% Customize the IRF plots
options.saveas_strng    = 'zerossigns-robust';
options.shocksnames     = {'ADshock','ASshock','MPshock'}; % 
% options.add_irfs        = bvar4.irzerosign_ols(indx_var,:,indx_sho,:);
plot_all_irfs_(irfs_to_plot,options)

return

%% IRFs with signs and higher-order moment (HOM) restriction
%--------------------------------------------------------------------------
% WARNING: THIS PART IS TIME  CONSUMING
% It takes about 22 minutes with a personal computer with an Intel(R)
% Core(TM) Ultra 7 155H 1.40 GHz processor and with 32.0 GB of RAM installed. 
%--------------------------------------------------------------------------
clear

load Data
% select the Euro area variables of interest 
y = [IPI CORE UNR NFCRATE Euribor1Y IT_DE_SPREAD];
optionsHOM.varnames = {'IP','Core','UNR','NFC Spread','Euribor 1y','It-De 5y'};

lags = 6;
optionsHOM.K   = 500;
optionsHOM.hor = 48;

% identify a spread shock in the EA:
% index of the shocks of interest
indx_sho              = 6;

% sign restrictions
optionsHOM.signs{1} = ['y(6,1:2,'  num2str(indx_sho) ')>0']; % IT-DE spread goes up
optionsHOM.signs{2} = ['y(4,1:2,'  num2str(indx_sho) ')>0']; % EBP goes up

% HOM restrictions, order of declaration of the setting: 
% 1. dimension: Shocks index, 
% 2. dimension: Estimator type: 
% sample mom = 1, 
% robust estimator = 2, 3, 4 (various types see skewness_.m and kurtosis_m
% for more details).
% 3. dimension: Moment order: skewness = 1; kurtosis = 2
% impose moderate skewness and kurtosis to the shock.
optionsHOM.hmoments{1} = ['hm('  num2str(indx_sho) ',4,2)> 0.5']; % kurtosis
optionsHOM.hmoments{2} = ['hm('  num2str(indx_sho) ',4,1)> 0.18']; % skewness - 
% Pearson's Coefficient of Skewness #2: 3*(mean-med)/s > 3*0.2
% A value between -1 and -0.5 or between 0.5 and 1 indicates a moderately
% skewed distribution. A value between -1.5 and -1 or between 1 and 1.5
% indicates a highly skewed distribution. A value less than -1.5 or greater
% than 1.5 indicates an extremely skewed distribution.    
optionsHOM.robust_bayes = 2;
bvarHOM       = bvar_(y,lags,optionsHOM);

optionsHOM.shocksnames   = {'Spread-KrtSkwn'};
optionsHOM.saveas_strng  = 'KrtSkwn';
optionsHOM.saveas_dir    = './irfs_plt';

indxf = find(isnan(squeeze(bvarHOM.irhmomsign_draws(1,5,indx_sho,:)))==0);
irfs_to_plot2           = bvarHOM.irhmomsign_draws(:,:,indx_sho,indxf);
optionsHOM. nplots = [2 3];
plot_irfs_(irfs_to_plot2,optionsHOM)
