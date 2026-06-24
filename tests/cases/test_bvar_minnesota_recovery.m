function result = test_bvar_minnesota_recovery()
%TEST_BVAR_MINNESOTA_RECOVERY Ground-truth recovery check for bvar_.
%
% Simulates data from a known VAR(2) and checks that bvar_, run with a
% Minnesota prior, recovers a posterior mean of Phi reasonably close to
% the true coefficients used to generate the data.
%
% This is a CORRECTNESS test, not a regression test: it does not depend
% on any stored "golden master" file, so it can catch bugs that are
% already present in the current code (e.g. a transposed coefficient
% matrix, a sign error, swapped equations), not just future changes that
% break previously-working behavior.
%
% Tolerance is intentionally loose (0.35, applied with a floor of 1 in
% assertClose -- see helpers/assertClose.m) because: (a) the Minnesota
% prior deliberately shrinks estimates away from the unrestricted OLS
% value, and (b) 300 observations and K=50 draws is a small/fast sample
% chosen for test speed, not estimation precision. The goal here is to
% catch gross errors, not to benchmark estimator efficiency.

result.name = 'bvar_minnesota_recovery';

ny   = 3;
lags = 2;
T    = 300;
[y, Phi_true, ~] = simulate_var(ny, lags, T, 42);

options              = struct();
options.priors.name  = 'Minnesota';
options.K            = 50;   % kept below 100 so bvar_ does not pop up a GUI waitbar
options.hor          = 4;
options.fhor         = 1;
options.noprint      = 1;

BVAR = bvar_(y, lags, options);

Phi_hat = mean(BVAR.Phi_draws, 3);

[passed, message] = assertClose(Phi_hat, Phi_true, 0.35, 'Phi posterior mean vs true Phi');

result.passed  = passed;
result.message = message;

end
