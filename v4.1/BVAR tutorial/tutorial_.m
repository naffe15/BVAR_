%% BVAR tutorial
% Author:   Filippo Ferroni
% Date:     27/02/2020

clear all
close all
clc

addpath ..\..\cmintools\
addpath ..\..\v4.1\

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

%% Extra part 1/: Forecast Error variance decomposition with Cholesky (bvar1)
% compute the contribution of MP to the H-step ahead forecast error  of
% observables

% consider the meand of the posterior distribution 
Phi   = mean(bvar1.Phi_draws,3);
Sigma = mean(bvar1.Sigma_draws,3);
% index of the shocks of interest (shock to gs1)
indx_sho              = [3];   

% 2 year ahead forecsat error
hh      = 8;
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
optnsplt.saveas_dir    = '.\sdcmp_plt';
% limit the plot to a specific time window 
optnsplt.Tlim          = [2006 2012];
plot_sdcmp_(yDecomp,bvar4,optnsplt)


% %==========================================================================
% % without initial condition
% 
% AAA0 = AAA;
% AAA0(:,:,end) = []; % without initial condition / deterministic component
% 
% pplotvar                = {'logip','logcpi','gs1'};
% optnsplt.plotvarnames   = {'IP','CPI','Interest Rate'};
% 
% 
% ex_names_ = { {'Shck 1'};...    Supply
%     {'Shck 2'};...              Demand
%     {'Shck 3'};...              MP    
%     };
% tag       = {'Supply';
%             'Demand';
%             'MP';
%             'Other Shocks'};
% 
% 
% optnsplt.tags = 'noy0';
% 
% plot_shcks_dcmp_(pplotvar,ex_names_,tag,AAA0,bvar3,optnsplt)

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
% savefigure_pdf('.\irfs_plt\IRFUN3D');

%% %%=========================================================================
%%% FAVAR %%%
%%=========================================================================
keyboard
clear all

% load favar data (Quarterly)
load DataFAVAR
% y2 slow moving variables (YoY growth rates)
transf = 2;         % standardize y2
% extrac the first 3 Principal Components (PC)
nfac   = 3;         % number of PC
[~,fhat,Lambda,~,STD] = pc_T(y2,nfac,transf);

% y1 interest rate (TBILL3M)

% compressed slow moving variables first (PC) and append TBILL3M,
y = [fhat y1];

% 1. restrictions on the compressed variables, recursive identification 
lags   = 2; 
fabvar = bvar(y,lags);

% Rescale the IRF from (PC+y1) back to (y2+y1)
% PC are ordered frist
order_pc = 1;
C_       = rescaleFAVAR(STD,Lambda,size(y1,2),order_pc);

% construct an IRF for each draw and shock of interest.
% shocks of interest: MP (3 factor + interest rate)
indx_sho     = nfac + 1;       
for k = 1: fabvar.ndraws % iterate on draws
    fabvar.irX_draws(:,:,1,k) =  C_ * fabvar.ir_draws(:,:,indx_sho,k);
end

% Indentify the variables of interest for plots (real GDP and CORE PCE)
[~,indx_var] = ismember ({'GDPC96','JCXFE'},varnames_y2);
% Real GDP and CORE PCE and TBILL3M
irfs_to_plot = fabvar.irX_draws( indx_var, :, 1, :);

% % Customize the IRF plot
% % variables names for the plots
options.shocksnames   = {'MP'};  
options.varnames      = {'GDP','CORE PCE'};  
plot_irfs_(irfs_to_plot,options)

%# sign restrictions on the uncompressed variables. 
% agregate supply: GDP (+) GDP deflator (-). 
% Assume that AD is the frist shock
[~,indx_var] = ismember ({'GDPC96','GDPCTPI'},varnames_y2);

signrestriction{1} = ['y(' num2str(indx_var(1)) ',2:3,1)>0;'];
signrestriction{2} = ['y(' num2str(indx_var(2)) ',2:3,1)<0;'];

for k = 1 : fabvar.ndraws % iterate on draws
      Phi       = fabvar.Phi_draws(:,:,k);
      Sigma     = fabvar.Sigma_draws(:,:,k);
      [ir,Omeg] = iresponse_sign(Phi,Sigma,fabvar.hor,signrestriction,C_);
      fabvar.irXsign_draws(:,:,:,k) = ir;
end

[~,indx_var] = ismember ({'GDPC96','GDPCTPI','JCXFE'},varnames_y2);
indx_sho     = 1;       % shocks of interest
irfs_to_plot = fabvar.irXsign_draws(indx_var ,:,indx_sho,:);

options.saveas_dir    = './irfs_plt';   % folder
options.saveas_strng  = 'FaVAR';        % names of the figure to save
options.shocksnames   = {'AS'};         % name of the shock
options.varnames      = {'GDP','GDP Defl','CORE PCE'};  
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)


%% %%=========================================================================
%%% PREDICTION %%%
%%=========================================================================
keyboard
% Housekeeping
clear all
close all
clc
% load the data
load Data
% select the variables of interest for the forecast exercise
yactual = [IPI HICP CORE Euribor1Y M3 EXRATE];
% stop at August 2014 for estimation
in_sample_end = find(T==2014 + 7/12);
y             = yactual(1:in_sample_end,:);
TT            = T(1:in_sample_end);

%% 1\ 
% Unconditional forecasts using flat priors
lags         = 6;
% one year forecasts
options.fhor = 12;
b.var(1) 	 = bvar(y,lags,options);

%% 2\ 
% Unconditional forecasts using Minnesota priors with default values

options.priors.name = 'Minnesota';
b.var(2)		    = bvar(y,lags,options);

%% 3\ 
% Unconditional forecasts using Minnesota priors with hyperparameter values
% maximizing the log data density 

options.max_minn_hyper  = 1;
options.index_est       = 1:4;
options.max_compute     = 3; % sims
options.lb              = [0 0 0 0]; % setting the lower bound
b.var(3)                = bvar(y,lags,options);

%% 4\ 
% Forecasts conditional on the path of the short run interest rate
% (decreasing path). Again we use Minnesota priors with hyperparameter
% values that maximize the log data density    
 
options.max_minn_hyper      = 0;
options.minn_prior_tau      = b.var(3).prior.minn_prior_tau;
options.minn_prior_decay    = b.var(3).prior.minn_prior_decay;
options.minn_prior_lambda   = b.var(3).prior.minn_prior_lambda;
options.minn_prior_mu       = b.var(3).prior.minn_prior_mu;
% select the variables that is conditioned (Interest rate #4)
options.endo_index          = 4;
% impose a lower bound trajectory for the short term interest rate
options.endo_path   = Euribor1Y(in_sample_end+1:end);
c.var(1)            = bvar(y,lags,options);

%% 5\ 
% Forecasts conditional on the path of the short run interest rate
% (decreasing path) using only monetary policy shocks identified via
% Cholesky decomposition constructed using Minnesota priors with
% hyperparameter values  maximizing the log data density    

options.exo_index   = options.endo_index;
c.var(2)            = bvar(y,lags,options);

%% Forecast Plot: plot forecasts against actual values o


%Storing the mean forecasts for the plot
for i = 1 :3
    tmp                  = mean(b.var(i).forecasts.no_shocks,3);
    b.var(i).frcsts_plot = [y; tmp];
end
tmp0 = mean(c.var(1).forecasts.conditional,3);
cfrcsts_plot = [y; tmp0];
tmp                 = mean(c.var(2).forecasts.conditional,3);
cfrcsts2_plot       = [y; tmp];


KAPPA = 12;
ordering_transf = [1, 1, 1, 0, 1, 0, 0];
col= strvcat('b','r','g');
tmp_str = 'frcsts_plt';
mkdir(tmp_str);
tzero = find(T == 2014);
varnames = {'IPI' 'HICP' 'CORE' 'Euribor1Y' 'M3' 'EXRATE'};

if tzero < KAPPA + 1
    tzero = KAPPA + 1;
end
for jj= 1:length(varnames)
    figure('Name',varnames{jj})
    if ordering_transf(jj)
        yplot       = 100*(yactual(KAPPA+1:end,jj)-yactual(1:end-KAPPA,jj));
        ycfplot     =  100*(cfrcsts_plot(KAPPA+1:end,jj) - cfrcsts_plot(1:end-KAPPA,jj));
        ycf2plot    =  100*(cfrcsts2_plot(KAPPA+1:end,jj) - cfrcsts2_plot(1:end-KAPPA,jj));
    else
        yplot = yactual(KAPPA+1:end,jj);
        %         yfplot = outputkf.Ynansfill(KAPPA+1:end,jj);
        ycfplot = cfrcsts_plot(KAPPA+1:end,jj);
        ycf2plot = cfrcsts2_plot(KAPPA+1:end,jj);
    end
    for ji = 1 : 3
        hold on
        if ordering_transf(jj)
            yfor(:,ji) = 100*(b.var(ji).frcsts_plot(KAPPA+1:end,jj)-b.var(ji).frcsts_plot(1:end-KAPPA,jj));
        else
            yfor(:,ji) = b.var(ji).frcsts_plot(KAPPA+1:end,jj);
        end
    end
    up = max ([max(yfor(tzero - KAPPA :end,:)), max(ycf2plot(tzero - KAPPA :end,:)), max(yplot(tzero - KAPPA :end,:))]);
    down = min ([min(yfor(tzero - KAPPA :end,:)), min(ycf2plot(tzero - KAPPA :end,:)), min(yplot(tzero - KAPPA :end,:))]);
    shade([T(end-KAPPA)],[T(end)],[.85 .85 .85],up,down)
    for ji = 1 :3
        plot(T(tzero:end),yfor(tzero - KAPPA :end,ji),'LineWidth',2,'Color',col(ji))
    end
    %     plot(T(tzero:end),yfplot(tzero - KAPPA :end),'k-.','LineWidth',2)
    axis tight
    plot(T(tzero:end),ycfplot(tzero - KAPPA :end),'y*-','LineWidth',2)
    axis tight
    plot(T(tzero:end),ycf2plot(tzero - KAPPA :end),'cd-','LineWidth',2)
    axis tight
    plot(T(tzero:end),yplot(tzero - KAPPA :end),'ko-.','LineWidth',2)
    axis tight
    
    if jj == 1
        legend({'Forecast Period','1) Flat Prior','2) Minnesota','3) Minnesota Opt','4) Cond Forecasts','5) Cond Forecasts only MP','Actual'},...
            'Location','Best','FontSize',16)
    end
    set(    gcf,'position' ,[50 50 900 650])
    STR_RECAP = [tmp_str '\multiple_frscst_' varnames{jj}];
    savefigure_pdf(STR_RECAP);
    saveas(gcf,STR_RECAP,'fig');
    saveas(gcf,STR_RECAP,'eps');
end
close all

%% Plot Opt Minnesota Forecast with bands.

% select the forecast to plot (Opt Minnesota)
frcsts                    = b.var(3).forecasts.with_shocks;               
% declare the directory where the plots are saved (in .\frcsts_plt) - default no save
options.saveas_dir        = 'frcsts_plt';                    
% appearence of subplots 
options.nplots            = [3 2];                           
% start of the forecast plot - default first date in-sample data
options.time_start        = 2013;
% Transformations
% 12 = Year over Year percentage change with monthly data for IPI and HICP 
options.order_transform  = [12 12 12 0 12 0]; 
% Titles for subplot
options.varnames         = ...
              {'Industrial Production','HICP','CORE','Euribor 1Y','M3','EXRATE'};
% multiple credible set - default .68
options.conf_sig_2       = 0.9;
% add actual data if possible
options.add_frcst        = yactual;

plot_frcst_(frcsts,y,TT,options)


%% %=========================================================================
%%% MIXED FREQUENCY VAR %%%
%%=========================================================================
keyboard
% Housekeeping
clear all
close all
clc
% load the mixed frequency data
load DataMF
% select the variables of interest for the mixed freq exercise exercise
y = [GDP IPI HICP CORE Euribor1Y UNRATE];
% specify the # of lags
lags = 6;              

% T aggregation: the quarterly variable  
% xq(t) = 1/3( xm(t) + xm(t-1) + xm(t-2)) at least two lags are needed
options.mf_varindex     = 1; 
options.K               = 1000;    % # draws
options.priors.name     = 'Minnesota';
% estimate the var model
bvarmf                  = bvar(y,lags,options);

% sort 
sorty = sort(bvarmf.yfill,3);

tmp_str = 'mfvar_plt';
mkdir(tmp_str);

% plot: levels
figure('Name','Monthly EA GDP')
plot(T,y(:,1),'Linewidth',2,'marker','o','color','r'); hold on; plot(T,sorty(:,1,0.5*bvarmf.ndraws),'b','Linewidth',2)
legend('Quarterly data','MXFREQ VAR flow (x^q_t = 1/3( x^m_t +  x^m_{t-1} +  x^m_{t-2}))','Interpolation','location','SouthOutside')
set(    gcf,'position' ,[50 50 900 650])
savefigure_pdf([tmp_str '\MonthlyGDP']);

% plot: growth rates
figure('Name','Monthly EA GDP growth')
Dy = 100*(y(4:end,1)-y(1:end-3,1));
Y = sorty(:,1,0.5*bvarmf.ndraws);
DY = 100*(Y(4:end,1)-Y(1:end-3,1));
% DiY = 100*(bvarmf.yinterpol(4:end,1)-bvarmf.yinterpol(1:end-3,1));
plot(T(4:end),Dy,'ro');
hold on; plot(T(4:end),DY,'b','Linewidth',2)
% hold on; plot(T(4:end),DiY,'k-.','Linewidth',1.5)
legend('Quarterly data','MXFREQ VAR flow (x^q_t = 1/3( x^m_t +  x^m_{t-1} +  x^m_{t-2}))','Interpolation','location','SouthOutside')
set(    gcf,'position' ,[50 50 900 650])
savefigure_pdf([tmp_str '\MonthlyGDPgrowth']);


