function [estimator,cov_Hessian,ME1,ME2,ME1_std] = probit(Y,X,method_flag)

% Purpose: 
% Estimate Probit model and the marginal effects
% -----------------------------------
% Model:
% Yi* = Xi * Beta + ui , where normalized ui ~ N(0,1)
% Yi* is unobservable. 
% If Yi* > 0, we observe Yi = 1; If Yi* <= 0, we observe Yi = 0
% -----------------------------------
% Algorithm: 
% Maximum likelihood with analytic gradients
% -----------------------------------
% Usage:
% Y = dependent variable (n * 1 vector)
% X = regressors (n * k matrix)
% method_flag = numeric gradients (1), analytic gradients (2, default), analytic Hessian (3)
% -----------------------------------
% Returns:
% estimator = estimator corresponding to the k regressors
% cov_Hessian = covariance matrix of the estimator
% ME1 = marginal effects (average data)
% ME2 = marginal effects (individual average)
% ME_std = standard error of ME1
% -----------------------------------
% it  prints  out
% estimate, st.er, t-stat,  me1,  st-me1, me2
% ----------------------------------------
% Notes: 
% Probit model is subject to normalization.
% The variance of disturbances is set to 1, and a constant is added to X.
% 
% Written by Hang Qian, Iowa State University
% Contact me:  matlabist@gmail.com



if nargin < 2;  error('Incomplete data');   end
if nargin < 3;  method_flag = 2;    end

try
    JustTest = normcdf(1);
catch
    disp(' ')
    disp('Oooopse, Matlab statistics toolbox is not installed.')
    disp('You may download a compatibility package on my website.')
    disp('http://www.public.iastate.edu/~hqi/toolkit')
    error('Program exits.')
end
 
[nrow_x,ncol_x] = size(X);
[nrow_y,ncol_y] = size(Y);
if nrow_x < ncol_x;    X = X';    ncol_x = nrow_x;end
if nrow_y < ncol_y;    Y = Y';    ncol_y = nrow_y;end
if ncol_x < ncol_y;    Y_temp = Y;    Y = X;    X = Y_temp;end

[nobs,nreg] = size(X);

Y_unique = unique(Y);
if length(Y_unique) ~= 2
    disp('Ooopse, dependent variable should be binary.');
    disp('For illustration purpose, Y are grouped to binary data by its mean value')
end
Y = (Y > mean(Y));


add_const = 0;
const_check = all(X==1);
if ~any(const_check)
    disp('A constant terms is added to X due to normalization. ')
    X = [ones(nobs,1),X];
    nreg = nreg + 1;
    const_check = [true,const_check];
    add_const = 1;
end


c_initial = (X'*X)\(X'*Y)*(rand-0.5)*10;
switch method_flag
    case 1
        options = optimset('LargeScale','off','MaxFunEvals',1000,'Display','off');        
    case 2
        options = optimset('LargeScale','off','GradObj','on','MaxFunEvals',1000,'Display','off','DerivativeCheck','off');        
    case 3
        options = optimset('LargeScale','on','GradObj','on','Hessian','on','MaxFunEvals',2000,'MaxIter',1000,'Display','off','DerivativeCheck','off'); 
end

try
    [estimator,log_like] = fminunc(@(c)ML_PROBIT(c,Y,X,method_flag),c_initial,options);
catch
    method_flag = 1;
    [estimator,log_like] = fminsearch(@(c)ML_PROBIT(c,Y,X,method_flag),c_initial);
end


q = 2*Y - 1;
lambda = q .* normpdf(q.*(X*estimator)) ./ normcdf(q.*(X*estimator));
Hessian = -transpose(lambda.*(lambda+X*estimator)*ones(1,nreg).*X)*X;
fisher_info = -Hessian;
cov_Hessian = inv(fisher_info);

% Gradient_indiv = lambda*ones(1,nreg).*X;
% cov_BHHH = inv(Gradient_indiv'*Gradient_indiv);

std_c = sqrt(diag(cov_Hessian));
t_stat = estimator./std_c;



X_bar = mean(X);
ME1 = normpdf(X_bar*estimator) * estimator;
ME2 = mean(normpdf(X*estimator)) * estimator;

Jocobian = normpdf(X_bar*estimator) * (eye(nreg)-X_bar*estimator*(estimator*X_bar));
ME_var = Jocobian * cov_Hessian * Jocobian';
ME1_std = sqrt(diag(ME_var));


ME1(const_check) = NaN;
ME2(const_check) = NaN;
ME1_std(const_check) = NaN;
eval([char([81 72 49 61]),'[87 114 105 116 116 101 110 32 98 121];'])
eval([char([81 72 50 61]),'[32 72 97 110 103 32 81 105 97 110];'])



result = cell(nreg+1,7);
result(1,:) = {' ','Estimator','SE','t-stat','ME(avg. data)','ME_std','ME(ind. avg.)'};           
for m=1:nreg
    result(m+1,1) = {['C(',num2str(m-add_const),')']};
    result(m+1,2:7) = {estimator(m),std_c(m),t_stat(m),ME1(m),ME1_std(m),ME2(m)};
end
disp(' ')
disp(result)



dummy_check = all((X==1) | (X==0)) & (~const_check);
for m = 1:nreg
    if dummy_check(m)
        disp(['C(',num2str(m-add_const),') corresponds to binary data.'])
        X_dummy1 = X;
        X_dummy0 = X;
        X_dummy1(:,m) = ones(nobs,1);
        X_dummy0(:,m) = zeros(nobs,1);        
        Y_margin_modify_data_bar = normcdf(mean(X_dummy1)*estimator) - normcdf(mean(X_dummy0)*estimator);
        Y_margin_modify_indi_bar = mean(normcdf(X_dummy1*estimator) - normcdf(X_dummy0*estimator));
        disp(['Adjusted marginal effects (average data) is ',num2str(Y_margin_modify_data_bar)])
        disp(['Adjusted marginal effects (individual average) is ',num2str(Y_margin_modify_indi_bar)])
        disp(' ')
    end
end

disp(['Log likelihood: ',num2str(-log_like)])
disp(' ')

%fwrite(1, char([QH1,QH2,10,13]))
end

%-----------------------------------
%  ML_PROBIT
%-----------------------------------
function [log_like,Gradient_c,Hessian_c] = ML_PROBIT(c,Y,X,method_flag)

q = 2*Y-1;
probit_F = normcdf(q.*(X*c));
log_like = sum(log(probit_F));
log_like = -log_like;


if method_flag >= 2    
    lambda = q.*normpdf(q.*(X*c))./probit_F;
    Gradient_c = lambda'*X;
    Gradient_c = -Gradient_c;
end

if method_flag == 3        
    [~,nreg] = size(X); 
    Hessian_c = transpose(lambda.*(lambda+X*c)*ones(1,nreg).*X)*X;     
    Hessian_c = -Hessian_c;
end

end

