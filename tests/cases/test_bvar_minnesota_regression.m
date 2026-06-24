function result = test_bvar_minnesota_regression()
%TEST_BVAR_MINNESOTA_REGRESSION Regression check for bvar_ (Minnesota prior).
%
% Compares bvar_'s output on a fixed synthetic dataset against a frozen
% fingerprint stored in reference/bvar_minnesota_fingerprint.mat. If this
% test fails, something changed bvar_'s numerical output relative to the
% saved baseline -- that might be an intended fix (re-run
% tests/generate_reference.m deliberately, after checking the change by
% hand) or it might be an unintended regression worth investigating.
%
% Because bvar_ resets its own RNG seed internally (rng(999) near the top
% of bvar_.m), repeated calls with identical inputs are expected to be
% bit-for-bit reproducible, not just "close" -- so this test uses a tight
% tolerance (1e-6) rather than the loose tolerance used in the recovery
% test.
%
% Requires reference/bvar_minnesota_fingerprint.mat to exist; run
% tests/generate_reference.m once locally to create it.

result.name = 'bvar_minnesota_regression';

here     = fileparts(mfilename('fullpath'));
ref_file = fullfile(here, '..', 'reference', 'bvar_minnesota_fingerprint.mat');

if ~exist(ref_file, 'file')
    result.passed  = false;
    result.message = ['No reference file found. Run tests/generate_reference.m ' ...
                       'once locally to create reference/bvar_minnesota_fingerprint.mat.'];
    return
end

loaded = load(ref_file);
ref    = loaded.fingerprint;

ny   = 3;
lags = 2;
T    = 300;
[y, ~, ~] = simulate_var(ny, lags, T, 7);   % must match generate_reference.m exactly

options              = struct();
options.priors.name  = 'Minnesota';
options.K            = 50;
options.hor          = 4;
options.fhor         = 1;
options.noprint      = 1;

BVAR = bvar_(y, lags, options);

ir_mean = mean(BVAR.ir_draws, 4);

[p1, m1] = assertClose(mean(BVAR.Phi_draws, 3),   ref.Phi_mean,   1e-6, 'Phi posterior mean');
[p2, m2] = assertClose(mean(BVAR.Sigma_draws, 3), ref.Sigma_mean, 1e-6, 'Sigma posterior mean');
[p3, m3] = assertClose(squeeze(ir_mean(:,1,:)),   ref.ir_h1,      1e-6, 'IRF at horizon 1');
[p4, m4] = assertClose(BVAR.logmlike,             ref.logmlike,   1e-8, 'log marginal likelihood');

result.passed  = p1 && p2 && p3 && p4;
result.message = sprintf('%s || %s || %s || %s', m1, m2, m3, m4);

end
