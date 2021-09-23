%% BVAR tutorial: Bayesian Dynamic Factor model 
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     09/14/2021

clear all;
close all; 
clc;

addpath ../../cmintools/
addpath ../../bvartools/

rng('default');
rng(999);

% generate factors/observables

% Setting the parameters
% AR1 factors
Phi    = [0.75 0.0;0.1 0.8];
% Cov factors innovation
Sigma  = [1 -0.2; -0.2 2];
Q      =  chol(Sigma,'lower');
% factor loadings
Lambda = [1 0;
          0.3 1;
          0.5 -1.5;
          -1.1 -0.5;
          1.5 1.5;            
      5*randn(45,size(Phi,1))];  
% persistence of idisyncratic terrors
rho    = 0; 
% standard deviation of idisyncratic errors
sig    = 1;

% preallocating memory
% sample length
T      = 301;
f = zeros(T,size(Phi,1));
y = zeros(T,size(Lambda,1));
e = zeros(T,size(Lambda,1));
i = zeros(T,size(Phi,1));

for t = 1 : T-1
    i(t+1,:) = randn(size(f,2),1)';
    f(t+1,:) = (Phi*f(t,:)' + Q*i(t+1,:)')';    
    e(t+1,:) = (rho*e(t,:)' + sig*randn(size(Lambda,1),1))';
    y(t+1,:) = (Lambda*f(t+1,:)' + e(t+1,:)')';
end
f(1,:) = [];
y(1,:) = [];

subplot(2,1,1)
plot(y)
subplot(2,1,2)
plot(f)

% true IRFs
hor = 24;
true    = iresponse(Phi,eye(size(f,2)),hor,Q);
Lam_    = repmat(Lambda,1,1,size(f,2));
if exist('pagemtimes','builtin') == 5
    ir_true = pagemtimes(Lam_,true);
else
    for ff = 1 : size(f,2)
        ir_true(:,:,ff) = Lambda * true(:,:,ff);
    end
end


%% principal components

nfac   = 2;
transf = 0;

[~,pc,~,egv,~] = pc_T(y,nfac,transf);

% scree plot
figure,
plot(egv(1:10),'b','Linewidth',2);
grid on
title('Eigenvalues')
% STR_RECAP = [ dirname '/screeplot'];
% saveas(gcf,STR_RECAP,'fig');
% savefigure_pdf([STR_RECAP '.pdf']);



%% static factor model

% 2 static factors
nfac = size(f,2);
lags = 0; 
% priors
options.priors.F.Lambda.mean = 0;
options.priors.F.Lambda.cov  = 6;
% estimation command
[BSFM] = bdfm_(y,lags,nfac,options);
% consider a subset of draws
index = 5000:20:size(BSFM.Phi_draws,3);
% plot estimated and true factors
figure,
for gg=1:nfac
    subplot(nfac,1,gg)
    % draws of factors
    plot([squeeze(BSFM.f_draws(:,gg,index))],'Color',[0.7 0.7 0.7])
    hold on
    % true factor
    f_mean = mean(BSFM.f_draws(:,gg,index),3);
    sign_ = sign(corr(f_mean,f(:,gg)));
    plot(sign_ * f(:,gg),'b','Linewidth',2)
    hold on
    plot(f_mean,'k','Linewidth',2)
end
% save plot
dirname  = 'factor_plt';
mkdir(dirname)
STR_RECAP = [ dirname '/sfm'];
saveas(gcf,STR_RECAP,'fig');
savefigure_pdf([STR_RECAP '.pdf']);


%% dynamic factor model
clear options
nfac   = size(f,2);
lags   = round(size(Phi,2)/nfac);

% Priors options
% factors priors
options.priors.F.Phi.mean    = Phi;
options.priors.F.Phi.cov     = 10 * eye(size(Phi,1));
options.priors.F.Sigma.scale = Sigma;
options.priors.F.Sigma.df    = 4;
options.priors.F.Lambda.mean = 0;
options.priors.F.Lambda.cov  = 6;
% else Jeffrey prior on PhiSigma: options.priors.name = 'Jeffrey';
% idyosincratic error priors
options.priors.G.Sigma.scale = sig;
options.priors.G.Sigma.df    = 4;
% IRF options
% IV identification
options.proxy    = i(lags+2:end,1);
% Sign identification
options.signs{1} = 'y(4,1:3,1)<0'; %
options.signs{2} = 'y(5,1:3,1)>0'; %
% estimation command
[BDFM] = bdfm_(y,lags,nfac,options);
% consider a subset of draws
index = 5000:20:size(BDFM.Phi_draws,3);
% plot estimated and true factors
figure,
for gg=1:nfac
    subplot(nfac,1,gg)
    plot([squeeze(BDFM.f_draws(:,gg,index))],'Color',[0.7 0.7 0.7])
    hold on
    f_mean = mean(BDFM.f_draws(:,gg,index),3);
    sign_ = sign(corr(f_mean,f(:,gg)));
    % true factor
    plot(sign_ * f(:,gg),'b','Linewidth',2)
    hold on
    plot(f_mean,'k','Linewidth',2)
end
STR_RECAP = [ dirname '/dfm'];
saveas(gcf,STR_RECAP,'fig');
savefigure_pdf([STR_RECAP '.pdf']);

%%

Phi_m = mean(BDFM.Phi_draws(:,:,index),3);
Sig_m = mean(BDFM.Sigma_draws(:,:,index),3);
Ptt  = lyapunov_symm(Phi_m, Sig_m);
f_m  = mean(BDFM.f_draws(:,:,index),3);
figure,
for gg=1:nfac
    subplot(nfac,1,gg)
    plot([f_m(:,gg)-2*Ptt(gg,gg) f_m(:,gg)+2*Ptt(gg,gg)],'Color',[0.7 0.7 0.7],'Linewidth',2)
    hold on
    plot(f(:,gg),'b','Linewidth',2)
    hold on
    plot(f_m(:,gg),'k','Linewidth',2)
end

index = 5000:20:size(BDFM.Phi_draws,3);
Phi0 = Phi';
jj = 0;
figure('Name','Phi')
for gg=1: nfac
    for hh =1 : nfac
        jj = jj+1;
        subplot(nfac,nfac,jj)
        histogram(squeeze(BDFM.Phi_draws(gg,hh,index)),40)
        hold on
        plot([Phi0(gg,hh) Phi0(gg,hh)],[0 50],'r','Linewidth',2)
        xlim([-1 1])
    end
end
%

figure('Name','Sigma')
Sigma0 = Q*Q';
jj = 0;
for gg=1:2
    for hh =1 :2
        jj = jj+1;
        subplot(nfac,nfac,jj)
        histogram(squeeze(BDFM.Sigma_draws(gg,hh,index)),40)
        hold on
        plot([Sigma0(gg,hh) Sigma0(gg,hh)],[0 50],'r','Linewidth',2)
    end
end
% 
% figure('Name','Lambda')
% Lambda0 = Lambda;
% jj = 0;
% for gg=1:10
%     for hh =1 :2
%         jj = jj+1;
%         subplot(5,4,jj)
%         histogram(squeeze(BDFM.lambda_draws(gg,hh,index)),40)
%         hold on
%         plot([Lambda0(gg,hh) Lambda0(gg,hh)],[0 30],'r','Linewidth',2)
%     end
% end
% 
% figure('Name','sigma_g')
% sig0 = sig*ones(size(Lambda,1),1);
% for jj=1:10
%     subplot(3,4,jj)
%     histogram(squeeze(BDFM.sigma_draws(jj,index)),40)
%     hold on
%     plot([sig0(jj) sig0(jj)],[0 30],'r','Linewidth',2)
% end

%% Plot IRFs

index_var          = [1 2 4 5 15 20];
index_sho          = 1;

% some options:
% add the true IRF
options.add_irfs   = ir_true(index_var,:,index_sho);
% additional 90% HPD set
options.conf_sig_2 = 0.9;   
% additional 90% HPD set
options.nplots = [2 3];   
% variables names for the plots
options.varnames      = {'Var1','Var2','Var4','Var5','Var15','Var20'};  
% name of the directory where the figure is saved
options.saveas_dir    = './factor_plt';
% name of the figure to save
options.saveas_strng  = 'sign';
% sign restricted IRF
irfs_to_plot       = BDFM.irsign_draws(index_var,:,index_sho,:);
plot_irfs_(irfs_to_plot,options)

% IV IRF
irfs_to_plot_iv     = BDFM.irproxy_draws(index_var,:,index_sho,:);
% name of the figure to save
options.saveas_strng  = 'iv';
plot_irfs_(irfs_to_plot_iv,options)






