% Example_1_classical.m Inference with a  Flat  Prior
% Author:   Filippo Ferroni
% Date:     27/02/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of  optimal  lag  lengh of   a VAR
% Computation of responses  and  vairance  decomposition  due  to  
% monetary policy  shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clc
addpath ../../cmintools/
addpath ../../bvartools/

 
% load monthly data from Getler  and  Karadi (2015) data set
load DataGK
y = [logip logcpi gs1 ebp];


% 1) search  for  optimal  lag  length
minlag  =5;  % minum number  of  lags
maxlag = 24; % maximum number of lags
opt.K  = 1;  % generate only 1 draw from the posterior
% 
for nlags=minlag:maxlag
  BVAR = bvar_(y, nlags, opt);
  disp(['Number of lags' num2str(nlags)])
  disp(BVAR.InfoCrit)  
end	
pause;

% 2) run the BVAR with  GK  lag lenght and Cholesky decomposition
lags        = 12;       % number  of  lags
options.hor = 48;       % number  of  horizons
bvar1       = bvar_(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho              = 3;   
% Change the order of the variables for the plot 
% 1. gs1; 2. logcpi; 3. logip; 4. ebp
indx_var              = [3, 2, 1, 4];
% IRF to PLOT
irfs_to_plot           = bvar1.ir_draws(indx_var,:,indx_sho,:);

% Customize the plot
% variables names for the plots
options.varnames      = {'1 year rate','CPI','IP','EBP'};  
% name of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% name of the figure to save
options.saveas_strng  = 'flat';
% name of the shock
options.shocksnames   = {'MP'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;   
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)
pause;

 
% 3) calculate variance  decomposition
% using the mean of the posterior distribution 
Phi   = mean(bvar1.Phi_draws,3);
Sigma = mean(bvar1.Sigma_draws,3);
% index of the shocks of interest (shock to gs1)
indx_sho              = 3;   

% 1/2/4 years ahead FEVD
cc=1;
FEVD=zeros(size(y,2), size(y,2),3);
for  hh       = [12 24 48]
    FEVD(:,:,cc)  = fevd(hh,Phi,Sigma);
    cc=cc+1;
end

disp('%=====================================================%')
disp('% Forecast Error Variance Decomposition               %')
disp('% Percentage of 1/2/4 years ahead volatility          %')
disp('% explained  by a  MP shock                           %')
disp('    logip     logcpi   gs1        ebp ')
disp(squeeze(FEVD(:,indx_sho,:))')
disp('%                                                     %')
disp('%=====================================================%')



 