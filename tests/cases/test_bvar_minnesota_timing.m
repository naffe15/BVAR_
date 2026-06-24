function result = test_bvar_minnesota_timing()
%TEST_BVAR_MINNESOTA_TIMING Coarse smoke test: bvar_ shouldn't blow up.
%
% This is NOT a precise performance benchmark. Wall-clock time depends
% heavily on the machine running it, so the pass/fail threshold below is
% deliberately generous -- meant to catch a gross regression (e.g.
% someone reintroducing an O(n^4)-ish computation in the draw loop), not
% to certify a specific speed. The measured time is reported in the
% message every time, pass or fail, so you can watch it over time even
% though the threshold itself is loose.
%
% For an actual before/after measurement of the kron-vs-direct
% optimization in bvar_'s Phi draw, see tests/benchmark_kron_vs_direct.m
% instead -- it reports a relative speedup ratio (hardware-independent),
% which is far more informative than this test's pass/fail outcome.
%
% NOTE: the threshold (60s) below was picked without being able to run
% this myself -- I have no MATLAB/Octave in the environment I wrote this
% in. After your first run, look at the actual reported time and tighten
% THRESHOLD_SEC if you want this to be a more sensitive regression check.

result.name = 'bvar_minnesota_timing';

ny   = 6;
lags = 4;
T    = 250;
[y, ~, ~] = simulate_var(ny, lags, T, 5);

options              = struct();
options.priors.name  = 'Minnesota';
options.K            = 99;   % kept below 100 so bvar_ does not pop up a GUI waitbar
options.hor          = 8;
options.fhor         = 1;
options.noprint      = 1;

THRESHOLD_SEC = 60;

t0 = tic;
BVAR = bvar_(y, lags, options); %#ok<NASGU> -- only the timing matters here
elapsed = toc(t0);

result.passed = elapsed <= THRESHOLD_SEC;
if result.passed
    result.message = sprintf(...
        'ny=%d, lags=%d, K=%d ran in %.2fs (threshold %.0fs, deliberately generous)', ...
        ny, lags, options.K, elapsed, THRESHOLD_SEC);
else
    result.message = sprintf(...
        'ny=%d, lags=%d, K=%d took %.2fs, exceeding the %.0fs threshold -- possible major performance regression', ...
        ny, lags, options.K, elapsed, THRESHOLD_SEC);
end

end
