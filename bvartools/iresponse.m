function [ir]=iresponse(alpha,Sigma,hor,Omega,unit,fast)
%IRESPONSE Impulse responses of a VAR for a given rotation of the shocks.
%
%   ir = IRESPONSE(alpha,Sigma,hor,Omega)
%   ir = IRESPONSE(alpha,Sigma,hor,Omega,unit)
%   ir = IRESPONSE(alpha,Sigma,hor,Omega,unit,fast)
%
%   Inputs
%     alpha   (ny*lags+nx) x ny  VAR coefficients, lag-1 block first; rows
%                                after ny*lags (constant, exogenous) are ignored
%     Sigma   ny x ny            covariance of the reduced-form residuals
%     hor     scalar             number of horizons (horizon 1 is the impact)
%     Omega   ny x m             orthonormal rotation; eye(ny) gives the
%                                recursive (Cholesky) responses
%     unit    any value          if passed (and not empty), the Cholesky factor
%                                is rescaled to a unit diagonal, so shock j
%                                moves variable j by 1 on impact (before Omega
%                                is applied). Omit it, or pass [], for
%                                one-standard-deviation shocks.
%     fast    0 or 1 (default 1) algorithm for the recursion (same output):
%               1: propagates the (ny*lags x m) block F^k*G*Q*Omega,
%                  cost O(hor*(ny*lags)^2*m)
%               0: original algorithm, propagates the full power F^k,
%                  cost O(hor*(ny*lags)^3)
%                With ny = 60, lags = 12 the fast version is about 10 times
%                faster (see examples/Stress test/stress_test_bvar.m).
%
%   Output
%     ir      ny x hor x m       response of variable i at horizon k to shock j
%
%   If Sigma is not positive definite, a factor based on the LDL
%   decomposition of pinv(Sigma) replaces the Cholesky factor.
%
%   See also IRESPONSE_SIGN, IRESPONSE_SET, IRESPONSE_LONGRUN, PLOT_IRFS_.

%   Filippo Ferroni, 6/1/2015.

if nargin < 6 || isempty(fast)
    fast = 1;
end

N           = size(Sigma,1);
% ir          = zeros(N,hor,N); % variables, horizon, shock
ir          = zeros(N,hor,size(Omega,2)); % variables, horizon, shock/exogenous variable
[m , n]     = size(alpha);
if rem(m, n)==0
    lags = m/n;
else
    lags = floor((m-1)/n);
end
try
    [Q]         = chol(Sigma,'lower');
catch
    iS               = pinv(Sigma);
    [L, D, P]        = ldl(iS);
    S_inv_upper_chol = sqrt(D) * L' * P';
    Q = S_inv_upper_chol';
end


% units
if nargin < 5 || isempty(unit)
    unit = eye(N);
else
    unit = inv(diag(diag(Q)));
end

% 1 standard deviation increase (if unit = eye(N))
% 100 basis point increase (else)
Q = Q*unit;

% companion form
F       = [alpha(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
G       = eye(N * lags, N);
% compute IRFs
if fast == 1
    % B_k = F^(k-1) * G * Q * Omega, updated as B_(k+1) = F * B_k; the
    % responses are its first N rows (G' * B_k)
    B = G * Q * Omega;
    for k=1:hor
        ir(1:N,k,:) = B(1:N,:);
        if k < hor
            B = F * B;
        end
    end
else
    Fk      = eye(N * lags);
    for k=1:hor
        PHI         = Fk * G * Q * Omega;
        ir(1:N,k,:) = G' * PHI;
        Fk          = F * Fk;
    end
end
