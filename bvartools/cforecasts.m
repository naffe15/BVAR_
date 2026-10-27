function [sims_with_endopath,EPSn] = cforecasts(endo_path,endo_index,forecast_data,Phi,Sigma)
%CFORECASTS Forecast conditional on the path of some variables, using all shocks.
%
%   [sims_with_endopath,EPSn] = CFORECASTS(endo_path,endo_index,forecast_data,Phi,Sigma)
%
%   Draws one path of the VAR that matches the given path of the
%   conditioning variables. All the (Cholesky) shocks can move; among the
%   shock sequences that reproduce the path, the minimum-norm part is
%   combined with a random draw of the remaining directions.
%
%   Inputs
%     endo_path      fhor x nc    path of the conditioning variables
%     endo_index     1 x nc       their positions in y (row vector)
%     forecast_data  struct with .initval (lags x ny) and .xdata (fhor x nx);
%                    see FORECASTS
%     Phi            (ny*lags+nx) x ny   VAR coefficients
%     Sigma          ny x ny             residual covariance
%
%   Outputs
%     sims_with_endopath  fhor x ny   conditional forecast (one draw)
%     EPSn                fhor x ny   orthogonal (Cholesky) shocks that
%                                     generate it
%
%   See also CFORECASTS2, FORECASTS, BVAR_.
%
%   References: Waggoner, D. F. and T. Zha (1999), "Conditional Forecasts
%   in Dynamic Multivariate Models", Review of Economics and Statistics
%   81(4), 639-651. Jarocinski, M. (2010), "Conditional Forecasts and
%   Uncertainty about Forecast Revisions in Vector Autoregressions",
%   Economics Letters 108(3), 257-259.

%   Filippo Ferroni, 6/1/2015.

[fhor,Ncondvar] = size(endo_path);
Nvar            = size(Sigma,2);

if Ncondvar ~= length(endo_index)
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
err = reshape(err',length(endo_index)*fhor,1);

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
% return the shocks as fhor x ny, as in cforecasts2 (old: ny x fhor)
EPSn = EPSn';
