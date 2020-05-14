%% BVAR tutorial: IRF
% Author:   Filippo Ferroni
% Date:     27/02/2020

clear all
close all
clc
addpath ../../cmintools/
addpath ../../v4.1/

%% %=========================================================================
%%% CAUSALITY %%%
%%=========================================================================

% load the data
load DataGK
y = [logip logcpi gs1 ebp];

%% 1/ Cholesky

% run the BVAR
lags        = 12;
options.hor = 48;
bvar1       = bvar(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho              = [3];   
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


%% 2/ Signs

% specify the restrictions
options.signs{1} = 'y(3,1:3,1)>0'; % 1Y rate up in period 1 to 3
options.signs{2} = 'y(2,1:3,1)<0'; % CPI down in period 1 to 3 
% run the BVAR
bvar2             = bvar(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho         = [1];  
% IRF ot plot (notice that we change the order of the variables for the plot)
irfs_to_plot     = bvar2.irsign_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% names of the figure to save
options.saveas_strng  = 'signs';
% name of the shock
options.shocksnames   = {'MPtightening'}; % 
plot_irfs_(irfs_to_plot,options)


%% 3/ Signs and Narrative

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
bvar3             = bvar(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho         = [1];  
% IRF ot plot (notice that we change the order of the variables for the plot)
irfs_to_plot     = bvar3.irnarrsign_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% names of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% names of the figure to save
options.saveas_strng  = 'signsnarrative';
% name of the shock
options.shocksnames   = {'MPtightening'}; % 
plot_irfs_(irfs_to_plot,options)


%% 4/ Zeros and Signs

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
bvar4             = bvar(y,lags,options);

% IRFs of Interest
indx_sho     = [1:3];
irfs_to_plot = bvar4.irzerosign_draws(indx_var,:,indx_sho,:);

% Customize the IRF plots
options.saveas_strng    = 'zerossigns';
options.shocksnames     = {'ADshck','ASshck','MPshck'}; % 
plot_irfs_(irfs_to_plot,options)

%% 5/ Long Run Technology

% Housekeeping: remove previous identification settings
options = rmfield(options,'zeros_signs');
options = rmfield(options,'saveas_strng');
options = rmfield(options,'shocksnames');

% define the dataset for the identification of LR shock
y = [diff(logip) logcpi(2:end) gs1(2:end) ebp(2:end)];

% run the BVAR
options.long_run_irf   = 1; 
% options.zeros_signs{1} = 'yl(1,1)=0;';
bvar5                  = bvar(y,lags,options);

% IRF of Interest
% shock index
indx_sho     =  [1];
% Define the order of the variables for the plot 
% 1. D(logip); 2. logcpi; 3.  gs1; 4. ebp;
indx_var     = [1:4];
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
 

%% 6/ Proxy or IV 

% define the dataset for the identification of MP shock with IV
y = [gs1 logip logcpi ebp];

% load the instruments
[numi,txti,rawi] = xlsread('factor_data.csv','factor_data');
% use the same instrument as GK
options.proxy  = numi(:,4);

% run the BVAR
bvar6        = bvar(y,lags,options);

% IRF of Interest
% shock index
indx_sho     =  [1];
% Keep the same order of variables as in the estimation
% 1. gs1; 2. logip; 2. logcpi; 4. ebp;
indx_var     = [1:4];
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
% % consider the meand of the posterior distribution 
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


%% Extra part 1/: Forecast Error variance decomposition with Cholesky (bvar1)
% compute the contribution of MP to the H-step ahead forecast error  of
% observables

% consider the meand of the posterior distribution 
Phi   = mean(bvar1.Phi_draws,3);
Sigma = mean(bvar1.Sigma_draws,3);
% index of the shocks of interest (shock to gs1)
indx_sho              = [3];   

% 2 year ahead FEVD
hh      = 24;
FEVD    = fevd(hh,Phi,Sigma);

disp('%=====================================================%')
disp('% Forecast Error Variance Decomposition of MP shock   %')
disp('% Percentage of volatility explained by MP shock      %')
disp('    logip     logcpi   gs1        ebp ')
disp(FEVD(:,indx_sho)')
disp('%                                                     %')
disp('%=====================================================%')



%% Extra part 2/: Historical Decomposition 
% VAR with zero-sign restrictions

% uses the zero-sign restrictions average rotation
opts_.Omega         =  mean(bvar4.Omegaz,3); 
% by default it uses mean over posterior draws
[yDecomp,ierror]  = histdecomp(bvar4,opts_); 

% yDecomp = historical decomposition
% time, variable, shocks and initial condition
% ierror = structural innovation

%plot_sdcmp_(yDecomp,bvar4)

% Declare the names of the variables in the order they appear in the VAR
bvar4.varnames      = {'IP','CPI','Interest Rate','EBP'};
% select the variables for the plot of the historical decomposition
optnsplt.plotvar_    = {'Interest Rate','EBP'};
% select the shocks combination to report
optnsplt.snames_ = { {'Shck1','Shck2'};...    Combine Supply and Demand
    {'Shck3'};...              MP    
    {'Shck4'} ...              Other shock not identified in the VAR
    };
% declare the name of the shocks
optnsplt.stag_       = {'Supply+Demand';
            'MP';
            'Other Shocks';
            'Initial Condition'};
% name of the file to save
optnsplt.save_strng    = 'y0';
% define the time for the plot
optnsplt.time          = T(1+lags:end);
% define the directory where the plot is saved 
optnsplt.saveas_dir    = './sdcmp_plt';
% limit the plot to a specific time window 
optnsplt.Tlim          = [2006 2012];
plot_sdcmp_(yDecomp,bvar4,optnsplt)


%% %% Extra part 3/: Minnesota Priors IRF 


clear options
y = [logip logcpi gs1 ebp];
lags                = 12;
options.hor         = 48;
options.presample      = 12;
options.prior.name     = 'minnesota';
options.minn_prior_tau = 0.5;
bvar6                  = bvar(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho              = [3];   
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
% add the Cholesi IRF with flat prior
options.add_irfs      = squeeze(median(bvar1.ir_draws(indx_var,:,indx_sho,:),4));
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)

%% tricks 1/ Rolling windows
% housekeeping
clear all
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
    rollbvar   = bvar(yQ(timespan(w,:),:),lags,options);
    % normalize for the size of the shock
    norm           = rollbvar.ir_ols(indx_sho,1,indx_sho);
    % Collect the MP OLS response in the window
    rollIRF(:,:,w) =  squeeze(rollbvar.ir_ols(:,:,indx_sho))/norm;
end

figure('name','UNR')
tt = T(timespan(:,end)');
surf(tt,[1:options.hor],squeeze(rollIRF(2,:,:)))
axis tight

