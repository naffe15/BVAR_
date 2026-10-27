function [sims_with_endopath,EPSn] = cforecasts2(endo_path,endo_index,exo_index,forecast_data,Phi,Sigma,Omega,EPSi,epslags_)
%CFORECASTS2 Forecast conditional on the path of some variables, using selected structural shocks.
%
%   [sims_with_endopath,EPSn] = CFORECASTS2(endo_path,endo_index,exo_index,forecast_data,Phi,Sigma)
%   [...] = CFORECASTS2(...,Omega,EPSi,epslags_)
%
%   Only the structural shocks in exo_index are chosen to hit the path;
%   the other shocks are drawn at random (or given in EPSi).
%
%   Inputs
%     endo_path      fhor x nc   path of the conditioning variables
%     endo_index     1 x nc      their positions in y (row vector)
%     exo_index      1 x nc      structural shocks used to hit the path,
%                                as many as conditioning variables
%     forecast_data  struct with .initval and .xdata; see FORECASTS
%     Phi            (ny*lags+nx) x ny   VAR coefficients
%     Sigma          ny x ny             residual covariance
%     Omega          ny x ny     rotation of chol(Sigma,'lower')
%                                (default eye: recursive shocks)
%     EPSi           fhor x ny   values of the other structural shocks
%                                (default: N(0,1) draws)
%     epslags_       0 (default), 1, 2 or 3: use the response k periods
%                    after the shock, instead of the impact response, to
%                    compute the conditioning shocks
%
%   Outputs
%     sims_with_endopath  fhor x ny   conditional forecast
%     EPSn                fhor x ny   structural shocks, with the
%                                     conditioning ones in exo_index
%
%   See also CFORECASTS, FORECASTS, BVAR_.
%
%   References: Waggoner and Zha (1999); Maih, J. (2010), "Conditional
%   Forecasts in DSGE Models", Norges Bank Working Paper 2010/07; Dynare
%   reference manual.

%   Filippo Ferroni, 6/1/2015.

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

if Ncondvar ~= length(endo_index)
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
    e(:,t) = inv(RR(endo_index,exo_index))*(endo_path(t,:) - y(endo_index))';
    y = y + (RR(:,exo_index)*e(:,t))';
    lags_data(1:end-1,:) = lags_data(2:end, :);
    lags_data(end,:) = y;
    sims_with_endopath(t, :) = y;
end

EPSn  = EPSi;
EPSn(:,exo_index) = e';


