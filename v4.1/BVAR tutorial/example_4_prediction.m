%% BVAR tutorial: FORECASTS
% Author:   Filippo Ferroni
% Date:     27/02/2020

clear all
close all
clc

addpath ../../cmintools/
addpath ../../v4.1/

%% %%=========================================================================
%%% PREDICTION %%%
%%=========================================================================
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
    STR_RECAP = [tmp_str '/multiple_frscst_' varnames{jj}];
    savefigure_pdf(STR_RECAP);
    saveas(gcf,STR_RECAP,'fig');
    saveas(gcf,STR_RECAP,'eps');
end
close all

%% Plot Opt Minnesota Forecast with bands.

% select the forecast to plot (Opt Minnesota)
frcsts                    = b.var(3).forecasts.with_shocks;               
% declare the directory where the plots are saved (in .\frcsts_plt) - default no save
options.saveas_dir        = '/frcsts_plt';                    
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

