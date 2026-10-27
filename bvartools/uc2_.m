function [ytrend,ycycle,out] = uc2_(y,lags,opts)
%UC2_ Unobserved-components model with an I(2) trend, estimated by maximum likelihood.
%
%   [ytrend,ycycle,out] = UC2_(y,lags)
%   [ytrend,ycycle,out] = UC2_(y,lags,opts)
%
%   Model:  y(t) = a(t) + b(t)
%           a(t) = a1 a(t-1) + ... + ap a(t-p) + ea(t)   cycle
%           b(t) = c(t-1) + b(t-1) + eb(t)               trend
%           c(t) = c(t-1) + ec(t)                        trend growth
%   The likelihood comes from KFILTERNAN. The shock variances are
%   optimised in logs, so they stay positive; the AR coefficients of the
%   cycle are optimised in levels and can be negative.
%
%   Inputs
%     y       T x 1    series
%     lags    scalar   order p of the cycle
%     opts    struct (optional)
%       .phi, .sigma     starting values (defaults 0.5^j and 0.5)
%       .index_est       parameters to estimate (default all)
%       .max_compute     optimiser, as in BVAR_MAX_HYPER (default 3)
%       .lb, .ub         bounds for max_compute = 2
%       .one_sided       1 returns filtered instead of smoothed components
%       .noplot, .noprint, .time, .varname   output controls
%
%   Outputs
%     ytrend, ycycle   T x 1   smoothed (or filtered) trend and cycle
%     out              KFILTERNAN output plus .phi and .sigma at the mode
%
%   See also UC_OPT_HYPERPARA, UC2SS, KFILTERNAN, UC1_.

%   Filippo Ferroni, 6/1/2020.

objective_function = 'uc_opt_hyperpara';

npara = lags+3;

% autoregressive parameters
hyperpara = ones(lags+3,1);
for elle = 1 : lags
    hyperpara(elle) = 0.5^elle;
end
% STD
hyperpara(lags+1:lags+3) = 0.5; 


max_compute = 3;
lb          = -1e10*ones(length(hyperpara),1);
ub          = 1e10*ones(length(hyperpara),1);
% x0          = log(hyperpara);   % old: AR coefficients in logs
x0          = uc_to_opt(hyperpara,1:npara,lags);
index_fixed = [];
index_est   = 1:npara;
hyperpara_fidex =  [];
one_sided   = 0; 
noplot  = 0;
noprint = 0;
time = 1:length(y);
varname ={'Var1'};

if nargin > 2
    if isfield(opts,'phi') ==1
        % set the initial values for phi
        hyperpara(1:lags) = opts.phi;
        % x0              = log(hyperpara);   % old: AR coefficients in logs
        x0              = uc_to_opt(hyperpara,1:npara,lags);
    end    
    if isfield(opts,'sigma') ==1
        % set the initial values for the VARIANCES of the US shocks
        hyperpara(1+lags : end) = opts.sigma;
        % x0              = log(hyperpara);   % old: AR coefficients in logs
        x0              = uc_to_opt(hyperpara,1:npara,lags);
    end    
    if isfield(opts,'index_est') ==1
        clear x0 index_fixed index_est hyperpara_fidex;
        index_est       = opts.index_est;
        index_fixed     = setdiff(1:npara,index_est);
        % x0              = log(hyperpara(index_est));   % old: AR coefficients in logs
        x0              = uc_to_opt(hyperpara(index_est),index_est,lags);
        % hyperpara_fidex = log(hyperpara(index_fixed));       % old: AR coefficients in logs
        hyperpara_fidex = uc_to_opt(hyperpara(index_fixed),index_fixed,lags);
    end    
    if isfield(opts,'max_compute') ==1
        max_compute    = opts.max_compute;
    end
    if isfield(opts,'lb') ==1
        % lb    = log(opts.lb);   % old: AR coefficients in logs
        lb    = uc_to_opt(opts.lb,index_est,lags);
        if length(lb) ~=  length(x0)
            error('Mismatch between the size of lower bounds and the param vector');
        end
    end
    if isfield(opts,'ub') ==1
        % ub    = log(opts.ub);   % old: AR coefficients in logs
        ub    = uc_to_opt(opts.ub,index_est,lags);
        if length(ub) ~=  length(x0)
            error('Mismatch between the size of upper bounds and the param vector');
        end
    end
    if isfield(opts,'objective_function') ==1
        objective_function = opts.objective_function;       
    end
    if isfield(opts,'one_sided') ==1
        one_sided = opts.one_sided;       
    end
    if isfield(opts,'noplot') ==1
        noplot = opts.noplot;
    end
    if isfield(opts,'time') ==1
        time = opts.time;
        if length(time) ~= length(y)
            error('Time must be of the same size as y')
        end
    end
    if isfield(opts,'varname') ==1
        varname = opts.varname;
    end

end
    
    
options.index_est       = index_est;
options.index_fixed     = index_fixed;
options.hyperpara_fidex =  hyperpara_fidex;

options.state_space_model = 2;
options.only_logL    = 1; 
options.initialCond  = 2; 

options.aZero        = zeros(lags+2+1,1);
options.aZero(end-2) = y(1);
% options.aZero(end-1) = y(1);
options.aZero(end)   = y(1);

options.pZero        = 10*eye(length(options.aZero));


switch max_compute
    case 1 % unconstraint
        % Set default optimization options for fminunc.
        optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
        [xh,fh,~,~,~,H] = ...
            fminunc(objective_function,x0,optim_options,y,lags,options);
     %=====================================================================     
     case 2 % constraint 
        % Set default optimization options for fmincon.
        optim_options = optimset('display','iter', 'LargeScale','off', 'MaxFunEvals',100000, 'TolFun',1e-8, 'TolX',1e-6);
        [xh,fh,~,~,~,~,H] = ...  
            fmincon(objective_function,x0,[],[],[],[],lb,ub,[],optim_options,y,lags,options);               
     %=====================================================================
    case 3 % sims 
        crit = 10e-5;
        nit  = 10e-4;
        [fh, xh, ~, H, ~, ~, ~] = csminwel(objective_function,x0,.1*eye(length(x0)),[],crit,nit,y,lags,options);
     %=====================================================================   
    case 7 % Matlab's simplex (Optimization toolbox needed).
%         optim_options = optimset('display','iter','MaxFunEvals',30000,'MaxIter',10000,'TolFun',1e-3,'TolX',1e-3,'OutputFcn',@outsavefun);
        optim_options = optimset('display','iter','MaxFunEvals',30000,'MaxIter',10000,'TolFun',1e-3,'TolX',1e-3);
        [xh,fh,~,~] = fminsearch(objective_function,x0,optim_options,y,lags,options);       
        H = zeros(length(x0));
end

postmode    = zeros(1,5);
% postmode(options.index_est)   = exp(xh);   % old: AR coefficients in logs
postmode(options.index_est)   = uc_from_opt(xh,options.index_est,lags);
% postmode(options.index_fixed) = exp(options.hyperpara_fidex);   % old: AR coefficients in logs
postmode(options.index_fixed) = uc_from_opt(options.hyperpara_fidex,options.index_fixed,lags);
if noprint == 0 % Print the output for Minnesoty prior with dummy

   
    % processing the output of the maximization
    log_dnsty   = -fh;
%     JJ          = jacob_bvar(xh);
%     HH          =  JJ * H * JJ';

    disp('=================================================================');
    disp('  ');
    disp('** Initial Hyperpara Values and Log Density **');
    disp('  ');
    str = ['rownam = {'];
    for ww=1:lags
        str = [str '''phi' num2str(ww) ''','];
    end
    str = [str '''siga'',''sigb'',''sigc'',''log density''};'];
    eval(str);
    minus_log_dnsty_0     = uc_opt_hyperpara(x0,y,lags,options);
    x                     = [hyperpara; -minus_log_dnsty_0];
    for jj =1: length(x)
        X = sprintf('%s = %0.5g',rownam{jj},x(jj));
        disp(X)
    end
    disp('  ');
    disp('** Posterior Mode: (Minimization of -Log Density)  **');
    disp('  ');
    x      = [postmode'; log_dnsty];
    for jj =1: length(x)
        X = sprintf('%s = %0.5g',rownam{jj},x(jj));
        disp(X)
    end
% else % other maximization e.g. conjugate MNIW
%     postmode    = exp(xh);
%     % processing the output of the maximization
%     log_dnsty   = -fh;
%     JJ          = jacob_bvar(xh);
%     HH          =  JJ * H * JJ';
end

phi   = (postmode(1:lags));
sigma = diag((postmode(1+lags:3+lags)));

options.only_logL = 0; 
[~ ,out]  = kfilternan(phi,sigma,y,options);
out.phi = phi;
out.sigma= sigma;

if one_sided
    ytrend = out.filteredSt(:,end-2); 
    ycycle = out.filteredSt(:,1);
else
    ytrend = out.smoothSt(:,end-2); 
    ycycle = out.smoothSt(:,1);
end    

if noplot==0
    %figure;
    subplot(2,1,1);
    hold on
    title(varname)
    %   plot(tid,tauhat, 'LineWidth',1,'Color','blue');
    plot(time,y,'r-','LineWidth',2);
    plot(time,ytrend, 'k--','LineWidth',2); axis tight;
    hold off
    box off; %xlim([tid(1)-1 tid(end)+1]);
    legend('Level','Trend','Location','NorthWest');
    subplot(2,1,2);
%     tmpy = repmat(y,1,2)-tauCI;
    hold on
    %    plotCI(tid,tmpy(:,1),tmpy(:,2));
    %    plot(tid, (y-tauhat), 'LineWidth',1,'Color','blue');
    %    plot(tid, zeros(T,1),'-k','LineWidth',1);
    plot(time,ycycle, 'b','LineWidth',2); axis  tight;
    %    plot(tid, zeros(T,1),'-k','LineWidth',1);
    hold off
    box off; %xlim([tid(1)-1 tid(end)+1]); ylim([-10 8])
    legend('Cycle','Location','NorthWest');
    set(gcf,'Position',[100 100 800 300])  %% 1 by 2 graphs
    pause
    close all;
end


end


% =========================================================================
function x = uc_to_opt(p,pos,lags)
% Parameters -> optimiser scale. pos are the positions of p in the full
% vector [phi(1:lags) sigma(1:3)]: AR coefficients stay in levels (they
% can be negative), variances are taken in logs.
x      = p(:);
k      = pos(:) > lags;
x(k)   = log(x(k));
end

function p = uc_from_opt(x,pos,lags)
% Inverse of uc_to_opt.
p      = x(:)';
k      = pos(:)' > lags;
p(k)   = exp(p(k));
end
