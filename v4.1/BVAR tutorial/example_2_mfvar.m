%% BVAR tutorial: MF
% Author:   Filippo Ferroni
% Date:     27/02/2020

clear all
close all
clc

addpath ../../cmintools/
addpath ../../v4.1/

%% %=========================================================================
%%% MIXED FREQUENCY VAR %%%
%%=========================================================================
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


