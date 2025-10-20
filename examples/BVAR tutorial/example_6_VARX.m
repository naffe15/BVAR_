% BVAR Tutorial: estimation of  VARX model
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     01/05/2020, revision 20/02/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate a  VAR for  the  UK with US and DE interest  rate  as exogenous
% variable. Bayesian estimation. 
% Plot responses to  US IR shocks. Do an  historical decomposition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear;
warning off;
addpath ../../cmintools/
addpath ../../bvartools/

% load the data
load DataPooling
% Time span:  1978m1 to 2012m8
cnames = {'uk','us','jp','de'};   % names of the countries
Nc = length(cnames);
% Variable names: IPI, CPI, 1Y GOVT YIELD (LTR), Policy Rate (STR)
vnames = {'ipi','cpi','ltr','str'};
Nv = length(vnames);

T       = size(time,1)-1;
 
% pick US rate (2,4) as  exogenous in a  VAR with UK variables(1,1:4)
lags        = 4 ;
options.hor = 24;
options.K   = 5000;
options.priors.name = 'Conjugate';
y = demean(100*diff(log ([ipi_uk ,cpi_uk , ltr_uk , str_uk]) )) ; 
% exogenous variables
z = demean(100*diff(log ([str_us str_de])));
% include one lag for the exogenous variable
options.controls = lagX(z,[0:1]);

% % signs/narrative restrictions
% options.signs{1}     = 'y(1,1:3,1)>0'; % 
% options.signs{2}     = 'y(2,1:3,1)<0'; % 
% options.narrative{1} = 'v([3],1)>0'; 
% % zero/sign restrictions
% % 1) ad = aggregate demand disturbance [sign restrictions]
% options.zeros_signs{1}     = 'y(1,1)=1;';
% options.zeros_signs{end+1} = 'y(2,1)=1;'; 
% options.zeros_signs{end+1} = 'y(3,1)=1;';
% % 2) as = aggregate supply shock [sign restrictions]
% options.zeros_signs{end+1} = 'y(1,2)=1;';
% options.zeros_signs{end+1} = 'y(2,2)=-1;';

% estimate the VARX
bvar1       = bvar_(y,lags,options); 

%% IRF to PLOT
indx_sho              = 1 : 2; 
indx_var              = [1, 2, 3, 4];
irfs_to_plot           = bvar1.irx_draws(indx_var,:,indx_sho,:);

% variables names for the plots
options.varnames      = {'UK IP','UK CPI', 'UK Long rate', 'UK Short rate'};  
% name of the shock
options.shocksnames   = {'US Short rate shock','DE Short rate shock'};

options.conf_sig_2    = 0.95;           % additional 95% HPD set
options.nplots        = [2 4];          % plot appeareance
options.saveas_strng  = 'VARX';         % name  of  the figure  to  save
options.saveas_dir    = './VARX_plt';   % folder where  it  is  stored
% the plotting command
plot_all_irfs_(irfs_to_plot,options);
pause;


%% shock decomposition in terms of domestic shocks and exogenous varibles
% yDecomp contains the decomosition of the observable variables in terms of
% (in this order)
% 1. Shocks (default identification - recursive)
% 2. Exogenous variables (if any)
% 3. purely deterministic component
[yDecomp,ierror]  = histdecomp(bvar1); 

% Declare the names of the variables in the order they appear in the VAR
bvar1.varnames      = options.varnames;
% select the variables for the plot of the historical decomposition
optnsplt.plotvar_   = options.varnames;
% select the shocks combination to report
optnsplt.snames_ = { {'Shck1','Shck2','Shck3','Shck4'};...    Combine Domestic shocks
    {'Shck5','Shck7'};...              US SRT at time (t) and (t-1)    
    {'Shck6','Shck8'} ...              DE SRT at time (t) and (t-1)    
    };
% declare the name of the shocks
optnsplt.stag_       = {'DomesticShocks';
            'US STR';
            'DE STR';
            'Initial Condition'};
% name of the file to save
optnsplt.save_strng    = 'decomp_exo_endo';
% define the time for the plot
% Time:         1978m1 to 2012m7
TT            = 1978 : 1/12 : 2012+6/12;
optnsplt.time = TT(1+lags:end);
optnsplt.Tlim = [2006 2012+6/12];
% define the directory where the plot is saved 
optnsplt.saveas_dir    = './VARX_plt';
% limit the plot to a specific time window 
% optnsplt.Tlim          = [2006 2012];
plot_sdcmp_(yDecomp,bvar1,optnsplt)
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extension with  2 exogenous variables lagged
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  WARNING : THIS  IS  SLOW

 z = demean(100*diff(log ([cpi_us str_us] )));
 options.controls = lagX(z,[0:1]);
 options.K   = 1000;
 % Sign restrictions
 options.signs{1} = 'y(4,1:3,4)>0';
 options.signs{2} = 'y(2,1:3,4)<0';

 % estimate the VARX
 bvar2 = bvar_(y,lags,options);
 
 indx_sho = 4;
 indx_var = [1, 2, 3, 4];
 irfs_to_plot = bvar2.irsign_draws(:,:,indx_sho,:);
 
 options.varnames = {'UK IP','UK CPI', 'UK Long rate', 'UK Short rate'};  
 options.shocksnames = {'MP_UK'};
 options.nplots= [2 2];
 plot_irfs_(irfs_to_plot,options);


%% Example for Exogenonus block with USA and Canada data about inflation, interest rates and unemployment rate
% Canada data is treated as exogenous in this analysis.
close all; clc; clear;

load DataEx.mat;

addpath ../../cmintools/
addpath ../../bvartools/

load DataEx.mat;
y = [CPI_US UNRATE_US INTRATE_US];
z = [CPI_CA UNRATE_CA INTRATE_CA];
lags = 4;

% Bvar 
options.exogenous_block = z;
bvarer = bvar_(y,lags,options);

% IRF graphs: all z (4:6) responses to y (1:3) shocks
[n, H, ~, ndraws] = size(bvarer.ir_draws);
horizons = 0:H-1;

figure;
for iz = 1:3        % z1,z2,z3  -> variabkes 4,5,6
    for iy = 1:3    % y1,y2,y3  -> shock 1,2,3
        
        response_draws = squeeze(bvarer.ir_draws(3+iz, :, iy, :)); % H × ndraws
        
        mean_irf  = mean(response_draws, 2);
        lower_irf = quantile(response_draws, 0.16, 2);
        upper_irf = quantile(response_draws, 0.84, 2);

        subplot(3,3,(iz-1)*3 + iy);
        hold on;

        fill([horizons, fliplr(horizons)], ...
             [lower_irf' fliplr(upper_irf')], ...
             [0.8 0.8 1], 'EdgeColor', 'none');

        plot(horizons, mean_irf, '-k', 'LineWidth', 1.5);
        yline(0,'--','Color',[0.5 0.5 0.5]);

        title(sprintf('z_%d \\leftarrow shock y_%d', iz, iy));
        if iz==3, xlabel('Horizon'); end
        if iy==1, ylabel('IRF'); end
        grid on;
    end
end

sgtitle('Impulse Responses of z to shocks in y');


figure;
for iz = 1:3        % z1,z2,z3  -> shocks 4,5,6
    for iy = 1:3    % y1,y2,y3  -> variables 1,2,3
        
        response_draws = squeeze(bvarer.ir_draws(iy, :, iz+3, :)); % H × ndraws
        
        mean_irf  = mean(response_draws, 2);
        lower_irf = quantile(response_draws, 0.16, 2);
        upper_irf = quantile(response_draws, 0.84, 2);

        subplot(3,3,(iy-1)*3 + iz);
        hold on;

        fill([horizons, fliplr(horizons)], ...
             [lower_irf' fliplr(upper_irf')], ...
             [0.8 0.8 1], 'EdgeColor', 'none');

        plot(horizons, mean_irf, '-k', 'LineWidth', 1.5);
        yline(0,'--','Color',[0.5 0.5 0.5]);

        title(sprintf('y_%d \\leftarrow shock z_%d', iy, iz));
        if iz==3, xlabel('Horizon'); end
        if iy==1, ylabel('IRF'); end
        grid on;
    end
end

sgtitle('Impulse Responses of y to shocks in z');


