% Example_6_VARX.m 
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     01/05/2020, revision 14/12/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate a  VAR for  the  UK with  US interest  rate  as exogenous
% variable. Estimation  is  Bayesian. Responses to  US IR shocks are
% plotted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc;

addpath ../../cmintools/
addpath ../../bvartools/

% load the data
load DataPooling
% Time:         1978m1 to 2012m8
cnames = {'uk','us','jp','de'};   % names of the countries
Nc = length(cnames);
% Varible names:    IPI, CPI, 1Y GOVT YIELD (LTR), Policy Rate (STR)
vnames = {'ipi','cpi','ltr','str'};
Nv = length(vnames);

T       = size(time,1)-1;
 
 
% pick US rate (2,4) as  exogenous in a  VAR with UK variables(1:1;4)
lags        = 4 ;
options.hor = 24;
options.K   = 1000;
options.priors.name = 'Conjugate';
y= demean(100*diff(log ([ipi_uk ,cpi_uk , ltr_uk , str_uk]) )) ; 
options.controls = demean(100*diff(log (str_us)));

bvar1       = bvar_(y,lags,options); 

% IRF to PLOT
indx_sho              = 1; 
indx_var              = [1, 2, 3, 4];
irfs_to_plot           = bvar1.irx_draws(indx_var,:,indx_sho,:);


% variables names for the plots
options.varnames      = {'UK IP','UK CPI', 'UK Long rate', 'UK Short rate'};  
% name of the shock
 options.shocksnames   = {'US Short rate shock'};

options.conf_sig_2    = 0.95;    % additional 95% HPD set
options.nplots = [1 4];          % plot appeareance
options.saveas_strng  = 'VARX';         % name  of  the figure  to  save
options.saveas_dir    = './VARX_plt';   % folder where  it  is  stored
% the plotting command
plot_all_irfs_(irfs_to_plot,options);
