function [y, Phi_true, Sigma_true] = simulate_var(ny, lags, T, seed)
%SIMULATE_VAR Generate data from a VAR(p) with known coefficients.
%
%   [y, Phi_true, Sigma_true] = simulate_var(ny, lags, T, seed)
%
% Used by the "recovery" style tests: an estimator run on data drawn from
% a known data-generating process should recover parameters close to the
% truth. Unlike a regression test, this does not depend on any stored
% reference file, so it also catches bugs that may already be present in
% the code under test, not just future regressions.
%
% Inputs:
%   ny    - number of variables
%   lags  - VAR lag order (p)
%   T     - number of observations to return (after burn-in)
%   seed  - RNG seed (default 42)
%
% Outputs:
%   y          - T x ny simulated data
%   Phi_true   - (ny*lags+1) x ny true coefficient matrix, in the SAME
%                row layout that bvartools uses internally (see YXB_.m
%                and BVAR.Phi_draws): rows stacked by lag, lag-1 block
%                first ... lag-p block last, intercept as the final row.
%                So Y = X*Phi_true + E with X = [y_{t-1} ... y_{t-p} 1].
%   Sigma_true - ny x ny true innovation covariance

if nargin < 4
    seed = 42;
end

if isOctave() == 0
    rng(seed, 'twister');
else
    randn('state', seed);
    rand('state', seed);
end

burn = 200;

% --- random, stability-checked lag coefficient blocks B{1}..B{lags} ---
% B{L}(i,j) = effect of the lag-L value of variable i on variable j,
% i.e. the same convention as the columns of YXB_'s regressor matrix.
B = cell(1, lags);
for L = 1:lags
    B{L} = (0.6/L) * (0.4*eye(ny) + 0.15*randn(ny, ny));
end

% Companion matrix (column-vector state form) just to check stability.
Comp = zeros(ny*lags, ny*lags);
for L = 1:lags
    Comp(1:ny, (L-1)*ny+1:L*ny) = B{L}';
end
if lags > 1
    Comp(ny+1:end, 1:ny*(lags-1)) = eye(ny*(lags-1));
end
rho = max(abs(eig(Comp)));
if rho >= 0.93
    shrink = 0.85/rho;
    for L = 1:lags
        B{L} = B{L} * shrink;
    end
end

const_true = 0.1*randn(1, ny);

% Well-conditioned, non-diagonal innovation covariance.
Lraw       = 0.3*randn(ny, ny);
Sigma_true = eye(ny) + Lraw*Lraw';
Cs         = chol(Sigma_true, 'lower');

Ttot        = T + burn + lags;
y           = zeros(Ttot, ny);
y(1:lags,:) = 0.5*randn(lags, ny);
for t = lags+1:Ttot
    yt = const_true;
    for L = 1:lags
        yt = yt + y(t-L,:)*B{L};
    end
    yt = yt + (Cs*randn(ny,1))';
    y(t,:) = yt;
end
y = y(burn+lags+1:end, :);

Phi_true = zeros(ny*lags+1, ny);
for L = 1:lags
    Phi_true((L-1)*ny+1:L*ny, :) = B{L};
end
Phi_true(end,:) = const_true;

end
