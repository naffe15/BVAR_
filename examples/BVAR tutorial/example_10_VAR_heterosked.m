%% BVAR tutorial: Minnesota prior and heteroskedasticity weights 
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     09/14/2020, revised 16/12/2020

close all; clc;

addpath ../../cmintools/
addpath ../../bvartools/


%% call the data

if exist('DataCovid.mat') == 2
    load DataCovid
else    
    url = 'https://fred.stlouisfed.org/';
    c = fred(url); 
    series    = {'PAYEMS','UNRATE','PCE','INDPRO','CPIAUCSL','PCEPILFE'};
    startdate = datenum('12/01/1988', 'mm/dd/yyyy');
    enddate   = datenum('07/01/2020', 'mm/dd/yyyy');
    time      = 1988+11/12 : 1/12 : 2020+6/12;
    
    for kk=1:size(series,2)
        tmp        = fetch(c,series{kk},startdate,enddate);
        if strcmp('UNRATE',series{kk})==1
            y(:,kk)    = tmp.Data(:,2);
            y0(:,kk)    = tmp.Data(:,2);
            y1(:,kk)  =(y0(:,kk));
        else
            y0(:,kk)    = tmp.Data(:,2);
            y1(:,kk)    = log(y0(:,kk));
            % rebase the log variable, transform to a index (2020m1=100)
            y(:,kk)= y1(:,kk) * 100 / y1(find(time==2019),kk);
        end
    end
    save DataCovid y y1 y0 time
end

% y, y0, y1 are  in DataCovid.

%% Forecast post COVID-19 (July 2020 - last insample data). Minnesota

lags                 = 13;
options.priors.name  = 'Minnesota';
options.fhor         = 24;
options.K            = 1000;
bvar0                = bvar_(y,lags,options);

options.nplots       = [2 3];                           
% start of the forecast plot - default first date in-sample data
options.time_start   = 2019;
% Titles for subplot
options.varnames = {'PAYMS','UNRATE','PCE','INDPRO','CPIAUCSL','PCEPILFE'};
% multiple credible set - default .68
options.conf_sig_2   = 0.9;
plot_frcst_(bvar0.forecasts.with_shocks,y,time,options)

pause;

% Forecast post COVID-19 (July 2020 - last insample data). 
% Minnesota + Heteroskedasticity weights

lags = 13;
tstar  = find(time==2020) + 2; % pick march 2020
% scale the variables by factor >1 in the periods that
% characterize the COVID-19 induced recession
st                   = ones(size(y,1),1);
st(tstar:tstar+2 ,:) = [10 10 10]; % March, April, May
st(1:lags)           = []; 

options.heterosked_weights = st;
options.priors.name  = 'Minnesota';
options.fhor         = 24;
bvar1                = bvar_(y,lags,options);

options.add_frcst = [y; mean(bvar0.forecasts.no_shocks,3)];
plot_frcst_(bvar1.forecasts.with_shocks,y,time,options)

pause;


%% this does not work (get  singular matrix  in  the  optimization
%% Maximize over Minnesota hyperparameters and heterosked weights

hyperpara(1)    = 3;		  % tau
hyperpara(2)    = 0.5;		  % decay
hyperpara(3)    = 1;		  % lambda
hyperpara(4)    = 1;		  % mu
hyperpara(5)    = 2;		  % omega
hyperpara(6)    = 2; % s0: scale march 2020
hyperpara(7)    = 2; % s1: scale april 2020
hyperpara(8)    = 2; % s2: scale may   2020
% setting the options
options.index_est          = [1 6:8]; % hyper-parameter over which maximize
options.objective_function = 'bvar_opt_heterosked';
options.tstar              = find(time==2020) + 2; %march 2020

[postmode,logmlike,HH] = bvar_max_hyper(hyperpara,y,lags,options);

heterosked_esse             = postmode(2:end); % s0, s1, s2   
esse                        = ones(size(y,1),1);
esse(options.tstar : options.tstar+2 ,:)   = heterosked_esse;
esse(1:lags)                = []; 

options.heterosked_weights = esse;

options.minn_prior_tau = postmode(1);
bvar2                 = bvar_(y,lags,options);

plot_frcst_(bvar2.forecasts.with_shocks,y,time,options)


