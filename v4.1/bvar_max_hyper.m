function [postmode,log_dnsty,HH] = bvar_max_hyper(hyperpara,y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'bvar_max_hyper' maximizes the marginal likelihood over the
% hyperparameters of for the Minnesota prior  

% Inputs:
% - hyperpara, Minnesota hyperpara over which maximize the marginal
% likelihood
% - y, data columns variables
% - lags, lag order of the VAR
% - options, see below for details

% Output: mode, Hessina and marginal likelihood at the mode

% Filippo Ferroni, 6/1/2017
% Revised, 3/21/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_compute = 3;
options_yes = 0;
lb          = -1e10*ones(length(hyperpara),1);
ub          = 1e10*ones(length(hyperpara),1);
objective_function = 'bvar_opt_hyperpara';

x0          = log(hyperpara);
index_fixed = [];
index_est   = 1:5;
hyperpara_fidex =  [];

if nargin > 3
    options_yes = 1;    
    if isfield(options,'index_est') ==1
        clear x0 index_fixed index_est hyperpara_fidex;
        index_est       = options.index_est;
        index_fixed     = setdiff(1:5,index_est);
        x0              = log(hyperpara(index_est));
        hyperpara_fidex = log(hyperpara(index_fixed));    
    end    
    if isfield(options,'max_compute') ==1
        max_compute    = options.max_compute;
    end
    if isfield(options,'lb') ==1
        lb    = options.lb;
        if length(lb) ~=  length(x0)
            error('Mismatch between the lower bounds and the param vector');
        end
    end
    if isfield(options,'ub') ==1
        ub    = options.ub;
        if length(ub) ~=  length(x0)
            error('Mismatch between the upper bounds and the param vector');
        end
    end
    if isfield(options,'objective_function') ==1
        objective_function = options.objective_function;
    end
end

options.index_est       = index_est;
options.index_fixed     = index_fixed;
options.hyperpara_fidex =  hyperpara_fidex;


switch max_compute
    case 1 % unconstraint
        % Set default optimization options for fminunc.
        optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
        if isfield(options,'optim_opt')
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end
        [xh,fh,exitflag,output,grad,H] = ...
            fminunc(objective_function,x0,optim_options,y,lags,options);
     %=====================================================================     
     case 2 % constraint 
        % Set default optimization options for fmincon.
        optim_options = optimset('display','iter', 'LargeScale','off', 'MaxFunEvals',100000, 'TolFun',1e-8, 'TolX',1e-6);
        if isfield(options,'optim_opt') && options_yes
            eval(['optim_options = optimset(optim_options,' options.optim_opt ');']);
        end
        [xh,fh,exitflag,output,lamdba,grad,H] = ...  
            fmincon(objective_function,x0,[],[],[],[],lb,ub,[],optim_options,y,lags,options);               
     %=====================================================================
    case 3 % sims 
        crit = 10e-5;
        nit  = 10e-4;
        [fh, xh, gh, H, itct, fcount, retcodeh] = csminwel(objective_function,x0,.1*eye(length(x0)),[],crit,nit,y,lags,options);
        
    case 7 % Matlab's simplex (Optimization toolbox needed).
%         optim_options = optimset('display','iter','MaxFunEvals',30000,'MaxIter',10000,'TolFun',1e-3,'TolX',1e-3,'OutputFcn',@outsavefun);
        optim_options = optimset('display','iter','MaxFunEvals',30000,'MaxIter',10000,'TolFun',1e-3,'TolX',1e-3);
        if isfield(options,'optim_opt') && options_yes
            eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
        end       
        [xh,fh,exitflag,others] = fminsearch(objective_function,x0,optim_options,y,lags,options);       
        H = zeros(length(x0));
end



postmode    = zeros(1,5);
postmode(options.index_est)   = exp(xh);
postmode(options.index_fixed) = exp(options.hyperpara_fidex);



% processing the output of the maximization
log_dnsty   = -fh;
JJ          = jacob_bvar(xh,options);
HH          =  JJ * H * JJ';


disp('=================================================================');
disp('  ');
disp('** Initial Hyperpara Values and Log Density **');
disp('  ');
rownam = {'tau','decay','lambda','mu','omega','log density'};
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

end

% if nargin< 4
%     x0 = log(hyperpara);
%
%     [fh, xh, gh, H, itct, fcount, retcodeh] = csminwel('bvar_opt_hyperpara',x0,.1*eye(length(x0)),[],10e-5,1000,y,lags);
%     postmode    = exp(xh);
% 	options.index_est = [1:length(xh)];
% else
%     options.index_fixed = setdiff(1:5,options.index_est);
%     x0 = log(hyperpara(options.index_est));
%     options.hyperpara_fidex = log(hyperpara(options.index_fixed));
%
%     [fh, xh, gh, H, itct, fcount, retcodeh] = csminwel('bvar_opt_hyperpara',x0,.1*eye(length(x0)),[],10e-5,1000,y,lags,options);
%     postmode    = zeros(1,5);
%     postmode(options.index_est)   = exp(xh);
%     postmode(options.index_fixed) = exp(options.hyperpara_fidex);
% end
% % processing the output of the maximization
% log_dnsty   = -fh;
% JJ          = jacob_bvar(xh,options);
% HH          =  JJ * H * JJ';
%
%
