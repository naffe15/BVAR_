function [sims_with_endopath,EPSn] = cforecasts2(endo_path,endo_index,exo_index,forecast_data,Phi,Sigma,Omega,EPSi,epslags_)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'cforecasts' performs the forecast conditional on a path for one or more
% endogenous variables using some of teh STRUCTURAL shocks of the VAR
% references: Waggoner and Zha (1999), Maih (2010), Dynare reference manual
% Inputs:
% - endo_path, path for the endogenous variables 
% - endo_index, order of the endogenous variables with a contional path
% - exo_index, order of the structural VAR shocks to use
% - forecast_data, last data
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - Omega, Rotation for identification

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

if nargin < 7
    Omega = eye(Nvar);
end
if nargin < 8
    % drawing the shocks
    % 1 standard deviation increase (if unit = eye(N))
    EPSi                = randn(fhor,Nvar);
    EPSi(:,exo_index)   = zeros(fhor,length(exo_index));
end
if nargin < 9
    epslags_ = 0;
end

if size(Ncondvar,2) ~= length(endo_index)
    error('something went wrong: the number of path and the number of vars are not consistent.')
end
if size(endo_index,1) > size(endo_index,2)
    error('the index needs to be a row vector');
end
if size(endo_index,2) ~= size(exo_index,2)
    error('the number of conditioning shocks needs to coincide with the number of conditioning endo');
end

[A,B,C,const,Sigma,lags]   = var2ss(Phi,Sigma);
sims_with_endopath = nan(fhor,Nvar);
EPSn               = nan(fhor,Nvar);

% % unconditional forecasts without zero
% [sims_no_shock,~]   = forecasts(forecast_data,Phi,Sigma,fhor,lags);
%
% % restrictions
% % err                 = endo_path_deviation - sims_no_shock(:,endo_path_index);
% % err = zeros(Nvar*fhor);
% err = (endo_path - sims_no_shock(:,endo_index))';
% % err = reshape(err',size(endo_index,1)*fhor,1);

[C]                 = chol(Sigma,'lower');
R                   = C*Omega;
EPS                 = (R*EPSi')';

if epslags_ == 1
    tmp = A*B*R;
    RR  = tmp(1:Nvar,1:Nvar);

elseif epslags_ == 2
    tmp = A^2*B*R;
    RR  = tmp(1:Nvar,1:Nvar);
    
elseif epslags_ == 3
    tmp = A^3*B*R;
    RR  = tmp(1:Nvar,1:Nvar);       
else
    RR  = R;
end

% With shocks but exo_index
lags_data = forecast_data.initval;
e         = zeros(length(exo_index),fhor);

for t = 1 : fhor
    X = [ reshape(flipdim(lags_data, 1)', 1, Nvar*lags) forecast_data.xdata(t, :) ];
    % unconditional forecasts with all shocks but 'exo_index'
    shock  = EPS(t,:);
    y = X * Phi + shock;
    e(:,t) = inv(RR(endo_index,exo_index))*(endo_path(t,:) - y(endo_index));
    y = y + (RR(:,exo_index)*e(:,t))';
    lags_data(1:end-1,:) = lags_data(2:end, :);
    lags_data(end,:) = y;
    sims_with_endopath(t, :) = y;
end

EPSn  = EPSi;
EPSn(:,exo_index) = e';


