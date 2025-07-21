%% BVAR tutorial: Mixed Frequency BVAR
% Author:   Filippo Ferroni and Fabio Canova
% Date:     27/02/2020, revised 21/07/2025

%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%%  estimation of  a  M-Q mixed  frequency  VAR
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------

close all; clc; clear ;

addpath ../../cmintools/
addpath ../../bvartools/
 
% load the mixed frequency data
load DataMF
y = [GDP IPI HICP CORE Euribor1Y UNRATE]; % variables  to  be  used
lags = 6;                                 % number  of  lags 

% T aggregation: the quarterly variable  
% xq(t) = 1/3( xm(t) + xm(t-1) + xm(t-2)) at least two lags are needed
options.mf_varindex     = 1; 
options.K               = 1000;   
options.priors.name     = 'Minnesota';
options.noprint         = 1;
bvarmf                  = bvar_(y,lags,options); % estimate the var model

sorty = sort(bvarmf.yfill,3);   % sort the smoothed  state


tmp_str = 'mfvar_plt';
mkdir(tmp_str);

% plot: levels
figure('Name','Monthly EA GDP')
plot(T,y(:,1),'Linewidth',2,'marker','o','color','r'); hold on; 
plot(T,sorty(:,1,0.5*bvarmf.ndraws),'b','Linewidth',2)
legend('Quarterly data','MXFREQ VAR flow (x^q_t = 1/3( x^m_t +  x^m_{t-1} +  x^m_{t-2}))','location','SouthOutside')
set(    gcf,'position' ,[50 50 900 650])
xlim([min(T), 2025 + 7/12]) 
if strcmp(version('-release'),'2022b') == 0
    savefigure_pdf([tmp_str '\MonthlyGDP']);
end

% plot: growth rates
figure('Name','Monthly EA GDP growth')
Dy = 100*(y(4:end,1)-y(1:end-3,1));
Y = sorty(:,1,0.5*bvarmf.ndraws);
DY = 100*(Y(4:end,1)-Y(1:end-3,1));
% DiY = 100*(bvarmf.yinterpol(4:end,1)-bvarmf.yinterpol(1:end-3,1));
plot(T(4:end),Dy,'ro');
hold on; plot(T(4:end),DY,'b','Linewidth',2)
% hold on; plot(T(4:end),DiY,'k-.','Linewidth',1.5)
legend('Quarterly data','MXFREQ VAR flow (x^q_t = 1/3( x^m_t +  x^m_{t-1} +  x^m_{t-2}))','location','SouthOutside')
set(    gcf,'position' ,[50 50 900 650])
xlim([min(T), 2025 + 7/12]) 
if strcmp(version('-release'),'2022b') == 0
    savefigure_pdf([tmp_str '\MonthlyGDPgrowth']);
end


%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------
%% Trailing nowcast example %%
%%-------------------------------------------------------------------------
%%-------------------------------------------------------------------------

close all; clc; clear ;
% load the mixed frequency data
load DataMF
y = [GDP IPI HICP CORE Euribor1Y UNRATE]; % variables  to  be  used
lags = 6;                                 % number  of  lags 

% Estimate the MFVAR model last full data arraw December 2024
LastDataPoint           = find(T==2025)-1; % december 2024
options.mf_varindex     = 1;
options.K               = 1000;
options.priors.name     = 'Minnesota';
options.noprint         = 1;
bvarmf                  = bvar_(y(1:LastDataPoint,:),lags,options);

% Construct the nowcast dataset
Tstart = 1;
YNowCast = nan(LastDataPoint -Tstart+1 +3, size(y,2));
YNowCast(1:LastDataPoint,:) = y(Tstart:LastDataPoint,:);
yNowCast = repmat(YNowCast,1,1,3*5);
jj = 0; tmp = YNowCast;
for mm = 1 : 3 % month (jan, feb, mar)
    for vv = 2 : 6 % variable (IPI HICP CORE Euribor1Y UNRATE)
        jj = jj + 1;
        % creating a dataset with increasing amount of information
        yNowCast(:,:,jj) = tmp;
        % adding one more datapoint
        yNowCast(LastDataPoint -Tstart+1 + mm,vv,jj) = y(LastDataPoint+mm,vv);
        % storing the trailing data
        tmp = yNowCast(:,:,jj);
    end
end

% Run the nowcast
NowOptions.fast_kf = 1;
NowOptions.noprint = 1;
out = nowcast_bvar(yNowCast, bvarmf, NowOptions);
% fast_kf=1
% Elapsed time is 26.411855 seconds.
% fast_kf=0
% Elapsed time is 353.161598 seconds.

% plot nowcast against actual

true = y(LastDataPoint + 3,1);
for dd = 1: size(yNowCast,3)
    GDP_now(:,dd) = squeeze(out.NowCast(LastDataPoint+3,:,dd));
end
GDP_dn = prctile(GDP_now,5,1);
GDP_up = prctile(GDP_now,95,1);
GDP_md = prctile(GDP_now,50,1);

tmp_str = 'mfvar_plt';
mkdir(tmp_str);
x = 1:dd;

for xx = 1 : dd
    figure(1),
    set(    gcf,'position' ,[50 50 900 650])
    subplot(2,2,1)
    plot(dd,true,'kd','MarkerSize',12,'LineWidth',2);
    hold on;
    plot(x(1:xx),GDP_md(1:xx),'Color',[0 0 0],'LineWidth',1.5);
    fill([x(1:xx) fliplr(x(1:xx))], [GDP_up(1:xx) fliplr(GDP_dn(1:xx))], [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    title('GDP - 2025Q1')
    xlim([1 dd+0.1])
    xlabel('data releases (#)');
    ylabel('log levels');
    grid on;
    ylim padded;

    subplot(2,2,2)
    plot(dd,100*(true-y(LastDataPoint+3-12,1)),'kd','MarkerSize',12,'LineWidth',2);
    hold on;
    plot(x(1:xx),100*(GDP_md(1:xx)-y(LastDataPoint+3-12,1)),'Color',[0 0 0],'LineWidth',1.5);
    fill([x(1:xx) fliplr(x(1:xx))], [100*(GDP_up(1:xx)-y(LastDataPoint+3-12,1)) fliplr(100*(GDP_dn(1:xx)-y(LastDataPoint+3-12,1)))], [0.8 1 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    title('GDP annual growth rates - 2025Q1')
    xlim([1 dd+0.1])
    xlabel('data releases (#)');
    ylabel('percent');
    grid on;
    ylim padded;

    subplot(2,2,3)
    plot(dd,100*(true-y(LastDataPoint+3-3,1)),'kd','MarkerSize',12,'LineWidth',2);
    hold on;
    plot(x(1:xx),100*(GDP_md(1:xx)-y(LastDataPoint+3-3,1)),'Color',[0 0 0],'LineWidth',1.5);
    fill([x(1:xx) fliplr(x(1:xx))], [100*(GDP_up(1:xx)-y(LastDataPoint+3-3,1)) fliplr(100*(GDP_dn(1:xx)-y(LastDataPoint+3-3,1)))], [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.1);
    title('GDP quarterly growth rates - 2025Q1')
    xlim([1 dd+0.1])
    xlabel('data releases (#)');
    ylabel('percent');
    grid on;
    legend('actual','Median','90% CI');
    ylim padded;

end
if strcmp(version('-release'),'2022b') == 0
    savefigure_pdf([tmp_str '\GDP_Now']);
end
