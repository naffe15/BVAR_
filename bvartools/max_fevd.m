function Qbar  = max_fevd(i, h, j, Phi, Sigma, Kappa)
%MAX_FEVD Rotation in which a shock explains most of a variable's FEVD.
%
%   Qbar = MAX_FEVD(i,h,j,Phi,Sigma)
%   Qbar = MAX_FEVD(i,h,j,Phi,Sigma,Kappa)
%
%   Draws Kappa random rotations and keeps, for each target, the one in
%   which shock j(s) has the largest share of the forecast error variance
%   of variable i(s) at horizon h.
%
%   Inputs
%     i       1 x s     target variables
%     h       scalar    horizon
%     j       1 x s     target shocks, one per variable in i
%     Phi     (ny*lags+nx) x ny   VAR coefficients
%     Sigma   ny x ny   covariance of the reduced-form residuals
%     Kappa   scalar    number of random rotations (default 1000)
%
%   Output
%     Qbar    ny x ny x s   chosen rotations; the columns are reordered so
%                           that the target shock j(s) is in column 1
%
%   See also FEVD, GENERATEQ.

if nargin < 6
    Kappa = 1000;
end
if length(i) ~= length(j)
    error('each shock must have one target')
end

N    = size(Sigma,1);
crit = nan(Kappa,length(j));
Q    = nan(N,N,Kappa);
for k = 1 : Kappa
    Q(:,:,k) = generateQ(N);                % generate an orthonormal matrix
    FEVD     = fevd(h,Phi,Sigma,Q(:,:,k));  % compute the FEVD
    % Calculate the contribution of shock j, to variable i forecast error volatility, at horizon h
    for j1 = 1 : length(j)
        crit(k,j1) = FEVD(i(j1), j(j1));
    end
end
% Pick the maximum
Qbar = nan(N,N,length(j));
for j1 = 1 : length(j)
    [~,index] = max(crit(:,j1));
    tmp_notindx = setdiff([1:N] , j(j1));
    Omega = Q(:,:,index);
    Qbar(:,:,j1) = Omega(:,[j(j1) tmp_notindx]);
end

end
