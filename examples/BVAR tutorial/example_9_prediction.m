%% BVAR tutorial: Out  of sample  predictions
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     27/02/2020, revised  20/02/2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Calculation of  unconditional  forecasts  with  different  priors
% 2) calculations  of conditional  forecasts using  different  assumptions
%%gda
% 3) structural scenarios (Gibbs sampler) with different identification schemes
%%gda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear;

%%gda
part = 3;    % parts to run: 1 = forecasts (Sections 1-5), 2 = structural scenarios (Sections 1-3 and 6), 3 = both
%%gda

addpath ../../cmintools/
addpath ../../bvartools/

% load the data
load Data
% select the variables of interest for the forecast exercise
yactual = [IPI HICP CORE Euribor1Y M3 EXRATE];
% stop at December 2024 for estimation - convention 0 = January, ...
in_sample_end = find(T==2023 + 11/12);
y             = yactual(1:in_sample_end,:);
TT            = T(1:in_sample_end);

%% 1)  Unconditional forecasts with flat prior
lags         = 6;
options.fhor = 12;             % one year forecasts
b.var(1) 	 = bvar_(y,lags,options);


%% 2)  Unconditional forecasts with default Minnesota prior
options.priors.name = 'Minnesota';
b.var(2)		    = bvar_(y,lags,options);

%% 3) Unconditional forecasts with optimal Minnesota  prior
% maximizing the log data density 
options.max_minn_hyper  = 1;
options.index_est       = 1:4;
options.max_compute     = 3;         % sims' optimization routine
options.lb              = [0 0 0 0]; % setting the lower bound
b.var(3)                = bvar_(y,lags,options);

%% 4 Forecasts conditional on the short run interest rate
%%gda
if part ~= 2
%%gda
% (decreasing path). Use  optimal Minnesota prior.    
 
options.max_minn_hyper      = 0;
options.minn_prior_tau      = b.var(3).prior.minn_prior_tau;
options.minn_prior_decay    = b.var(3).prior.minn_prior_decay;
options.minn_prior_lambda   = b.var(3).prior.minn_prior_lambda;
options.minn_prior_mu       = b.var(3).prior.minn_prior_mu;
% select the variables that is conditioned (Interest rate #4)
options.endo_index          = 4;
% impose a lower bound trajectory for the short term interest rate
options.endo_path   = Euribor1Y(in_sample_end+1:end);
% % impose a lower bound trajectory for the short term interest rate and M3
% options.endo_index  = [4 5];
% options.endo_path   = [Euribor1Y(in_sample_end+1:end) M3(in_sample_end+1:end)];
c.var(1)            = bvar_(y,lags,options);

%% 5) Forecasts conditional on the short run interest rate
% (decreasing path) using only monetary policy shocks identified via
% Cholesky decomposition and optimal Minnesota prior.
%  
options.exo_index   = options.endo_index;
c.var(2)            = bvar_(y,lags,options);

%% Forecast Plot: plot forecasts against actual values 

% Store the mean forecasts for the plot
for i = 1 :3
    tmp                  = mean(b.var(i).forecasts.no_shocks,3);
    b.var(i).frcsts_plot = [y; tmp];
end
tmp0 = mean(c.var(1).forecasts.conditional,3);
cfrcsts_plot = [y; tmp0];
tmp                 = mean(c.var(2).forecasts.conditional,3);
cfrcsts2_plot       = [y; tmp];


KAPPA = 12;
ordering_transf = [0, 1, 1, 0, 1, 0, 0];
%col= strvcat('b','r','g');
col= char('b','r','g');
tmp_str = 'frcsts_plt';
mkdir(tmp_str);
tzero = find(T == 2023 + 6/12);
varnames = {'IPI' 'HICP' 'CORE HICP' 'Euribor1Y' 'M3' 'EXRATE'};

if tzero < KAPPA + 1
    tzero = KAPPA + 1;
end
yfor=zeros(size(CORE,1)-KAPPA,3);
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
    shade(T(end-KAPPA),T(end),[.85 .85 .85],up,down)
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
    saveas(gcf,STR_RECAP,'fig');
    if strcmp(version('-release'),'2022b') == 0
        savefigure_pdf(STR_RECAP);
        saveas(gcf,STR_RECAP,'eps');
    end
end
pause;
close all

%% Plot Optimal Minnesota Forecast with bands.

% select the forecast to plot (Minnesota)
frcsts                    = b.var(3).forecasts.with_shocks;               
% declare the directory where the plots are saved (in .\frcsts_plt) - default no save
%%gda
% relative path: with the leading slash of the original ('/frcsts_plt') mkdir fails in the root folder
options.saveas_dir        = 'frcsts_plt';
%%gda
% appearence of subplots
options.nplots            = [3 2];                           
% start of the forecast plot - default first date in-sample data
options.time_start        = 2023 + 6/12;
% Transformations
% 12 = Year over Year percentage change with monthly data for IPI and HICP 
options.order_transform  = [0 12 12 0 12 0]; 
% Titles for subplot
options.varnames         = ...
              {'Industrial Production','HICP','CORE','Euribor 1Y','M3','EXRATE'};
% multiple credible set - default .68
options.conf_sig_2       = 0.9;
% add actual data if possible
options.add_frcst        = yactual;

plot_frcst_(frcsts,y,TT,options)


%%gda
end    % if part ~= 2
%% 6) Structural scenarios with different identification schemes
% Targets: realized paths of 2024, as in Sections 4-5. Each scenario is run with
% soft and hard uncertainty. Spare shocks are added to the driving shocks.

if part ~= 1
fhor = 12;
IPI_ = 1; HICP_ = 2; CORE_ = 3; RATE_ = 4; M3_ = 5;     % position in y
tgt  = in_sample_end + (1:fhor);
tol  = std(diff(y));                                    % tolerance of the soft paths

options_scn                   = struct();
options_scn.fhor              = fhor;
options_scn.priors.name       = 'Minnesota';
options_scn.minn_prior_tau    = b.var(3).prior.minn_prior_tau;
options_scn.minn_prior_decay  = b.var(3).prior.minn_prior_decay;
options_scn.minn_prior_lambda = b.var(3).prior.minn_prior_lambda;
options_scn.minn_prior_mu     = b.var(3).prior.minn_prior_mu;

% 1) Euribor1Y path, Cholesky
scenarios(1).name           = 'Euribor1Y path, Cholesky';
scenarios(1).id_opts        = struct();
scenarios(1).paths          = struct('var',RATE_,'periods',1:fhor,'value',yactual(tgt,RATE_));
scenarios(1).driving_shocks = [RATE_ IPI_ HICP_];
scenarios(1).uncertainty    = 'preserve_unconditional';    % tolerance = unconditional variance of the path

% 2) Euribor1Y path and M3 at months 6 and 12, zero and sign restrictions
% monetary policy shock (4): no impact on IPI, HICP, CORE, Euribor1Y up, M3 down; money shock (5): M3 up, Euribor1Y up
scenarios(2).name           = 'Euribor1Y path and M3 at two dates, zeros and signs';
scenarios(2).id_opts.zeros_signs = {'ys(1,4)=0;','ys(2,4)=0;','ys(3,4)=0;','y(4,4)=1;','y(5,4)=-1;', ...
                                    'y(5,5)=1;','y(4,5)=1;'};
scenarios(2).paths          = [struct('var',RATE_,'periods',1:fhor,'value',yactual(tgt,RATE_)), ...
                               struct('var',M3_,'periods',[6 12],'value',yactual(tgt([6 12]),M3_))];
scenarios(2).driving_shocks = [RATE_ M3_ IPI_ HICP_ CORE_];
scenarios(2).uncertainty    = 'user';
scenarios(2).Omega_f        = diag([repmat(tol(RATE_)^2,fhor,1); repmat(tol(M3_)^2,2,1)]);

% 3) HICP path, sign restrictions
% demand shock (1): IPI, HICP, CORE up; supply shock (2): IPI up, HICP, CORE down
% HICP only: with HICP and CORE together the chain can become explosive
scenarios(3).name           = 'HICP path, signs';
scenarios(3).id_opts.signs  = {'y(1,1:3,1)>0','y(2,1:3,1)>0','y(3,1:3,1)>0', ...
                               'y(1,1:3,2)>0','y(2,1:3,2)<0','y(3,1:3,2)<0'};
scenarios(3).paths          = struct('var',HICP_,'periods',1:fhor,'value',yactual(tgt,HICP_));
scenarios(3).driving_shocks = [1 2 RATE_ M3_];
scenarios(3).uncertainty    = 'user';
scenarios(3).Omega_f        = 4*tol(HICP_)^2*eye(fhor);    % inflated variance (standard deviation 2 times tol, about the unconditional one at month 12)

% hard twins
for k = 1:3
    scenarios(3+k)             = scenarios(k);
    scenarios(3+k).name        = [scenarios(k).name ' (hard)'];
    scenarios(3+k).uncertainty = 'hard';
    scenarios(3+k).Omega_f     = [];
end

%% Run the scenarios
scn = cell(numel(scenarios),1);
for i = 1:numel(scenarios)
    opts_run = options_scn;
    fn = fieldnames(scenarios(i).id_opts);
    for k = 1:numel(fn)
        opts_run.(fn{k}) = scenarios(i).id_opts.(fn{k});
    end
    opts_run.scenario.name           = scenarios(i).name;
    opts_run.scenario.horizon        = fhor;
    opts_run.scenario.paths          = scenarios(i).paths;
    opts_run.scenario.driving_shocks = scenarios(i).driving_shocks;
    opts_run.scenario.uncertainty    = scenarios(i).uncertainty;
    opts_run.scenario.Omega_f        = scenarios(i).Omega_f;
    scn{i} = bvar_(y,lags,opts_run);
end

%% Plot the scenarios
% red line = target path, dotted gray = unconditional forecast;
% plausibility q: 0.5 = same as the unconditional forecast, 1 = implausible
popt           = struct();
popt.varnames  = {'IPI','HICP','CORE','Euribor1Y','M3','EXRATE'};
popt.add_frcst = [y; median(b.var(3).forecasts.no_shocks,3)];
plot_frcst_(scn,y,TT,popt);
end    % if part ~= 1
%%gda
%% 
