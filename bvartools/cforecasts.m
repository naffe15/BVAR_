function [sims_with_endopath,EPSn] = cforecasts(endo_path,endo_index,forecast_data,Phi,Sigma)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'cforecasts' performs the forecast conditional on a path for one or more
% endogenous variables using all the shocks of the VAR
% references: Waggoner and Zha (1999) and Jarocinski (2010)

% Inputs:
% - endo_path, path for the endogenous variables 
% - endo_index, order of the endogenous variables with a contional path
% - forecast_data, last data
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks

% Output:
% - sims_with_endopath
% 1st dimension:   horizon 
% 2nd dimension:   variable 
% - EPSn, shocks generating the path

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[fhor,Ncondvar] = size(endo_path);
Nvar            = size(Sigma,2);

if size(Ncondvar,2) ~= length(endo_index)
    error('something went wrong: the number of path and the number of vars are not consistent.')
end
if size(endo_index,1) > size(endo_index,2)
    error('the index needs to be a row vector');
end
Nres = size(endo_index,2)*fhor;

[F,~,~,const,~,lags]     = var2ss(Phi,Sigma);

% endo_path_deviation = endo_path - repmat(const',1,fhor);

% no shocks forecasts
[sims_no_shock,~]   = forecasts(forecast_data,Phi,Sigma,fhor,lags);

% restrictions
% err                 = endo_path_deviation - sims_no_shock(:,endo_path_index);
% err = zeros(Nvar*fhor);
err = endo_path - sims_no_shock(:,endo_index);
err = reshape(err',size(endo_index,1)*fhor,1);

% 1 standard deviation increase (if unit = eye(N))
[C] = chol(Sigma,'lower');

R    = zeros(Nvar * fhor);
tmp0 = [];
for ff = 1 :fhor
    tmp  = F^(ff-1);
    tmp0 = [tmp(1:Nvar,1:Nvar) * C tmp0];
    R((ff-1)*Nvar + 1 : (ff)*Nvar, 1: size(tmp0,2) ) = tmp0;
    index(ff,:) = (ff-1)*Nvar + endo_index;
end
index = reshape(index',Nres,1);
% index= [4 10];
Rtilde = R(index,:);

[U,D,V] = svd(Rtilde);
V1 = V(:,1:Nres);
V2 = V(:,Nres+1:end);

eps = V1*inv(D(1:Nres,1:Nres))*U'*err + V2 * randn(size(R,1)-Nres,1);

% orthogonal schoks
EPSn = reshape(eps,Nvar,fhor);
% reduced form shocks
EPS  = C*EPSn;

[~,sims_with_endopath] = forecasts(forecast_data,Phi,Sigma,fhor,lags,EPS');
