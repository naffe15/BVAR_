% Example_1_classical.m Inference with a  Flat  Prior
% Author:   Filippo Ferroni
% Date:     27/02/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Estimation of  optimal  lag  lengh of   a VAR
% 2) Computation of responses  and  variance  decomposition  due  to  
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
  disp(['Number of lags ' num2str(nlags)])
  disp(BVAR.InfoCrit)  
end	
pause;

% 2) run the BVAR with  GK  lag lenght and Cholesky decomposition
lags        = 12;       % number  of  lags
options.hor = 48;       % number  of  horizons
bvar1       = bvar_(y,lags,options);

% Define the IRF of Interest

% Change the order of the variables for the plot 
% 1. gs1; 2. logcpi; 3. logip; 4. ebp
indx_var              = [3, 2, 1, 4];
% index of the shocks of interest; orthogonalized var
indx_sho              = indx_var;   
% IRF to PLOT
irfs_to_plot           = bvar1.ir_draws(indx_var,:,indx_sho,:);

% Customize the plot
% variables names for the plots
options.varnames      = {'1 year rate','CPI','IP','EBP'};  
% name of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% name of the figure to save
options.saveas_strng  = 'flat_allshocks';
% name of the shock
options.shocksnames   = {'eR','eCPI','eIP','eEBP'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;   
% normalize to a unitary change in the shock (default is STD options.normz = 0)
options.normz = 1;
% the plotting command
plot_all_irfs_(irfs_to_plot,options)
pause;

 
% 3) calculate variance  decomposition
% using the mean/median of the posterior distribution 
Phi   = mean(bvar1.Phi_draws,3);
Sigma = mean(bvar1.Sigma_draws,3);
Phi1   = median(bvar1.Phi_draws,3);
Sigma1 = median(bvar1.Sigma_draws,3);

% index of the shocks of interest (shock to gs1)
indx_sho              = 3;   

% 1/2/4 years ahead FEVD ( using  the mean  or  the  median)
cc=1;
FEVD=zeros(size(y,2), size(y,2),3); FEVD1=FEVD;
for  hh       = [12 24 48]
    FEVD(:,:,cc)  = fevd(hh,Phi,Sigma);  % mean
    FEVD1(:,:,cc)  = fevd(hh,Phi1,Sigma1);  % median
    cc=cc+1;
end

disp('%=====================================================%')
disp('% Forecast Error Variance Decomposition               %')
disp('% Percentage of 1/2/4 years ahead volatility          %')
disp('% Explained  by a  Monetary Policy shock              %')
disp('    logip     logcpi    gs1        ebp ')
disp(' Mean')
disp(squeeze(FEVD(:,indx_sho,:))')
disp('Median')
disp(squeeze(FEVD1(:,indx_sho,:))')
disp('%                                                     %')
disp('%=====================================================%')

% computing the distribution of the FEVD
options.conf_sig =0.95;
bvar1 = fevd_(12,bvar1,options);
disp('%=====================================================%')
disp('% Forecast Error Variance Decomposition               %')
disp('% Percentage of one year ahead volatility             %')
disp('% Explained  by a  Monetary Policy shock              %')
disp('    logip     logcpi    gs1        ebp ')
disp('                                                       ')
disp('Upper HPD')
disp(squeeze(bvar1.fevd.up(:,indx_sho))')
disp('Median')
disp(squeeze(bvar1.fevd.median(:,indx_sho))')
disp('Lower HPD')
disp(squeeze(bvar1.fevd.low(:,indx_sho))')
disp('%                                                     %')
disp('%=====================================================%')

