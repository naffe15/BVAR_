function [postmode,log_dnsty,HH] = bvar_max_hyper(hyperpara,y,lags,options)
%BVAR_MAX_HYPER Maximises the marginal likelihood over prior hyperparameters.
%
%   [postmode,log_dnsty,HH] = BVAR_MAX_HYPER(hyperpara,y,lags)
%   [postmode,log_dnsty,HH] = BVAR_MAX_HYPER(hyperpara,y,lags,options)
%
%   Optimises over the log of the hyperparameters, so they stay positive,
%   and prints the starting values and the mode.
%
%   Inputs
%     hyperpara  1 x n   starting values; for the default objective
%                        [tau decay lambda mu omega]
%     y          T x ny  data
%     lags       scalar  number of lags
%     options    struct (optional)
%       .index_est           positions of hyperpara to estimate (default
%                            all); the others stay fixed
%       .max_compute         1 fminunc, 2 fmincon (uses .lb, .ub),
%                            3 csminwel (default), 7 fminsearch
%       .lb, .ub             bounds, in levels, for max_compute = 2
%       .objective_function  name of the objective (default
%                            'bvar_opt_hyperpara'; also 'bvar_opt_plr',
%                            'bvar_opt_heterosked')
%       .hyperpara_names     names used in the printout
%       plus the options of BVAR_ML
%
%   Outputs
%     postmode   mode of the hyperparameters; with the default objective
%                all five, otherwise only the estimated ones
%     log_dnsty  log marginal likelihood (objective) at the mode
%     HH         Hessian of minus the log density at the mode, with respect
%                to the estimated hyperparameters in levels, for every
%                optimiser (zeros for fminsearch). inv(HH) approximates
%                their covariance.
%
%   See also BVAR_OPT_HYPERPARA, BVAR_OPT_PLR, BVAR_OPT_HETEROSKED, BVAR_ML,
%   CSMINWEL.

%   Filippo Ferroni, 6/1/2017.

max_compute = 3;
nhpara      = length(hyperpara);
lb          = -1e10*ones(nhpara,1);
ub          = 1e10*ones(nhpara,1);
dummy_      = 1;
x0          = log(hyperpara);
index_fixed = [];
index_est   = 1:nhpara;
hyperpara_fidex =  [];
mh_         = 0;
hyperpara_names = {};

objective_function = 'bvar_opt_hyperpara';

if nargin > 3
    if isfield(options,'index_est') ==1
        clear x0 index_fixed index_est hyperpara_fidex;
        index_est       = options.index_est;
        index_fixed     = setdiff(1:nhpara,index_est);
        x0              = log(hyperpara(index_est));
        hyperpara_fidex = log(hyperpara(index_fixed));
    end
    if isfield(options,'max_compute') ==1
        max_compute    = options.max_compute;
    end
    if isfield(options,'lb') ==1
        lb    = log(options.lb);
        if length(lb) ~=  length(x0)
            error('Mismatch between the size of lower bounds and the param vector');
        end
    end
    if isfield(options,'ub') ==1
        ub    = log(options.ub);
        if length(ub) ~=  length(x0)
            error('Mismatch between the size of upper bounds and the param vector');
        end
    end
    if isfield(options,'objective_function') ==1
        objective_function = options.objective_function;
        dummy_ =0;
    end
    if isfield(options,'hyperpara_names') ==1
        hyperpara_names = options.hyperpara_names;
        if length(hyperpara_names) ~=  nhpara
            error('Mismatch between the size of hyperpara_names and the hyperpara vector');
        end
    end
end


options.index_est       = index_est;
options.index_fixed     = index_fixed;
options.hyperpara_fidex =  hyperpara_fidex;


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
% HH below is the Hessian of minus the log density in levels, whatever the
% optimiser: csminwel returns an inverse Hessian, so invert it first
if max_compute == 3
    H = H \ eye(size(H));
end

if dummy_ % Print the output for Minnesoty prior with dummy
    postmode    = zeros(1,5);
    postmode(options.index_est)   = exp(xh);
    postmode(options.index_fixed) = exp(options.hyperpara_fidex);
    
    % processing the output of the maximization
    log_dnsty   = -fh;
    JJ          = jacob_bvar(xh);
    % HH          =  JJ * H * JJ';   % old: covariance rule, wrong for a Hessian
    iJJ         = JJ \ eye(size(JJ));
    HH          = iJJ * H * iJJ';    % Hessian in levels: inv(J)' * H * inv(J)
    
    disp('=================================================================');
    disp('  ');
    disp('** Initial Hyperpara Values and Log Density **');
    disp('  ');

    rownam = {'tau','decay','lambda','mu','omega'};              
    for pp = 1 : nhpara-5                                 
        rownam{end+1} = sprintf('pandemic phi(%d)', pp);           
    end
    rownam{end+1} = 'log density'; 
    
    minus_log_dnsty_0     = bvar_opt_hyperpara(x0,y,lags,options);
    x               = [hyperpara'; -minus_log_dnsty_0];
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
else % other maximization e.g. conjugate MNIW
    postmode    = exp(xh);
    % processing the output of the maximization
    log_dnsty   = -fh;
    JJ          = jacob_bvar(xh);
    % HH          =  JJ * H * JJ';   % old: covariance rule, wrong for a Hessian
    iJJ         = JJ \ eye(size(JJ));
    HH          = iJJ * H * iJJ';    % Hessian in levels: inv(J)' * H * inv(J)
    if isempty(hyperpara_names) == 0 % Print the output when the hyperparameters are named
        fullmode                    = hyperpara;
        fullmode(options.index_est) = postmode;

        disp('=================================================================');
        disp('  ');
        disp('** Initial Hyperpara Values and Log Objective **');
        disp('  ');
        rownam = [hyperpara_names(:)', {'log objective'}];
        minus_log_dnsty_0 = feval(objective_function,x0,y,lags,options);
        x = [hyperpara(:); -minus_log_dnsty_0];
        for jj =1: length(x)
            X = sprintf('%s = %0.5g',rownam{jj},x(jj));
            disp(X)
        end
        disp('  ');
        disp('** Mode: (Minimization of -Log Objective)  **');
        disp('  ');
        x = [fullmode(:); log_dnsty];
        for jj =1: length(x)
            X = sprintf('%s = %0.5g',rownam{jj},x(jj));
            disp(X)
        end
    end
end


end