function result = test_bvar_penalized()
%TEST_BVAR_PENALIZED Smoke test for the Ridge/Lasso/ElasticNet path in bvar_.
%
% This exercises the 'for pp = 1:3' block in bvar_.m that was rewritten to
% remove eval()-based dynamic field assignment (BVAR.(name).field = ...
% instead of eval(['BVAR.' name '.field = ...'])). It checks that the
% rewritten code still runs and produces a BVAR.<name> struct with the
% expected fields, finite values, and the right shapes -- a SMOKE test,
% not a precision check on the penalized estimators themselves.
%
% Ridge needs no extra toolbox and is always checked. Lasso/ElasticNet
% need MATLAB's Statistics and Machine Learning Toolbox (the 'lasso'
% function); if it isn't on the path, those two are SKIPPED -- reported
% in the message, but not counted as a failure, since that's a fact about
% the environment, not a bug in bvar_. (bvar_ itself already raises a
% clear error in this case: 'Cannot estimate VAR with Lasso: matlab stat
% toolbox needed' -- this test doesn't need to re-check that behavior.)

result.name = 'bvar_penalized_smoke';

ny   = 3;
lags = 2;
T    = 200;
[y, Phi_true, ~] = simulate_var(ny, lags, T, 11);

messages   = {};
all_passed = true;

% --- Ridge: no extra toolbox needed ---
options              = struct();
options.priors.name  = 'Minnesota';
options.K            = 50;
options.noprint      = 1;
options.Ridge        = struct('lambda', 0.05);

BVAR = bvar_(y, lags, options);

[ok, msg]      = check_penalized_block(BVAR, 'Ridge', ny, lags, Phi_true);
all_passed     = all_passed && ok;
messages{end+1} = msg;

% --- Lasso / ElasticNet: only if the Statistics Toolbox's lasso() exists ---
has_lasso = (exist('lasso') == 2); %#ok<EXIST> -- matches the check bvar_.m itself uses

if has_lasso
    options2             = struct();
    options2.priors.name = 'Minnesota';
    options2.K           = 50;
    options2.noprint     = 1;
    options2.Lasso       = struct('lambda', 0.05);
    BVAR2 = bvar_(y, lags, options2);
    [ok2, msg2]     = check_penalized_block(BVAR2, 'Lasso', ny, lags, Phi_true);
    all_passed      = all_passed && ok2;
    messages{end+1} = msg2;

    options3              = struct();
    options3.priors.name  = 'Minnesota';
    options3.K            = 50;
    options3.noprint      = 1;
    options3.ElasticNet   = struct('lambda', 0.05, 'alpha', 0.5);
    BVAR3 = bvar_(y, lags, options3);
    [ok3, msg3]     = check_penalized_block(BVAR3, 'ElasticNet', ny, lags, Phi_true);
    all_passed      = all_passed && ok3;
    messages{end+1} = msg3;
else
    messages{end+1} = 'Lasso/ElasticNet: SKIPPED (Statistics and Machine Learning Toolbox not found)';
end

result.passed  = all_passed;
result.message = join_(messages, ' | ');

end

% --------------------------------------------------------------------------
function [ok, msg] = check_penalized_block(BVAR, name, ny, lags, Phi_true)
% Checks BVAR.(name) exists with the expected fields, sizes, and finite
% values. Also reports (for information only -- this is not a pass/fail
% criterion) the max absolute coefficient error vs. the data's true Phi,
% since penalized estimators are expected to be biased/shrunk and a
% loose recovery check isn't meaningful here the way it is for bvar_
% Minnesota -- see test_bvar_minnesota_recovery.m for that style of test.

if ~isfield(BVAR, name)
    ok  = false;
    msg = sprintf('%s: FAILED (BVAR.%s missing)', name, name);
    return
end

blk             = BVAR.(name);
required_fields = {'Phi', 'e', 'Sigma', 'ir', 'InfoCrit'};
present         = isfield(blk, required_fields);

if ~all(present)
    ok  = false;
    msg = sprintf('%s: FAILED (missing fields: %s)', name, join_(required_fields(~present), ', '));
    return
end

size_ok = isequal(size(blk.Phi), [ny*lags+1, ny]) ...
    && all(isfinite(blk.Phi(:))) ...
    && all(isfinite(blk.Sigma(:))) ...
    && all(isfinite(blk.ir(:)));

if ~size_ok
    ok  = false;
    msg = sprintf('%s: FAILED (wrong size or non-finite values in Phi/Sigma/ir)', name);
    return
end

err = max(abs(blk.Phi(:) - Phi_true(:)));
ok  = true;
msg = sprintf('%s: OK (fields present, finite, max abs coef error vs true Phi = %.3g, expected to be biased)', name, err);

end

% --------------------------------------------------------------------------
function s = join_(parts, delim)
% Minimal strjoin replacement, kept local for MATLAB/Octave portability.
if isempty(parts)
    s = '';
    return
end
s = parts{1};
for i = 2:numel(parts)
    s = [s delim parts{i}];
end
end
