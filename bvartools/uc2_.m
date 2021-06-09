function [ytrend,ycycle,out] = uc2_(y,lags,opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from the UC of the form 
% a(t)  = a1 a(t-1) + ... + ap a(t-p) + ea(t )    [transition 1]
% b(t)  = c(t-1)    + b(t-1) + eb(t);             [transition 2]
% c(t)  = c(t-1)             + ec(t);             [transition 3]
% y(t)  = a(t)      + b(t);                       [transition 4]

% Filippo Ferroni, 6/1/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
x0          = log(hyperpara);
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
        x0              = log(hyperpara);
    end    
    if isfield(opts,'sigma') ==1
        % set the initial values for the VARIANCES of the US shocks
        hyperpara(1+lags : end) = opts.sigma;
        x0              = log(hyperpara);
    end    
    if isfield(opts,'index_est') ==1
        clear x0 index_fixed index_est hyperpara_fidex;
        index_est       = opts.index_est;
        index_fixed     = setdiff(1:npara,index_est);
        x0              = log(hyperpara(index_est));
        hyperpara_fidex = log(hyperpara(index_fixed));    
    end    
    if isfield(opts,'max_compute') ==1
        max_compute    = opts.max_compute;
    end
    if isfield(opts,'lb') ==1
        lb    = log(opts.lb);
        if length(lb) ~=  length(x0)
            error('Mismatch between the size of lower bounds and the param vector');
        end
    end
    if isfield(opts,'ub') ==1
        ub    = log(opts.ub);
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
postmode(options.index_est)   = exp(xh);
postmode(options.index_fixed) = exp(options.hyperpara_fidex);
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



