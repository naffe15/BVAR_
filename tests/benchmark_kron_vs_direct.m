%BENCHMARK_KRON_VS_DIRECT Quantify the bvar_ Phi-draw optimization.
%
% bvar_.m used to draw Phi from its matrix-normal conditional posterior by
% materializing kron(Sigma_lower_chol, XXi_lower_chol) -- a (ny*nk)x(ny*nk)
% matrix -- every single draw. It now uses the mathematically equivalent
% XXi_lower_chol*Z*Sigma_lower_chol' (via vec(A*Z*B') = kron(B,A)*vec(Z)),
% which never builds that matrix. The two give numerically identical draws
% (see tests/cases/test_bvar_minnesota_regression.m for the end-to-end
% check, which passed with a ~1e-17 relative error after the change).
%
% This script isolates just that one step and times it directly, across a
% few representative problem sizes, since the difference is hard to see
% from a full bvar_ run (which also spends time on OLS, IRFs, forecasts,
% etc. -- see the note printed at the end).
%
% Run manually:
%   run('<path to BVAR_>/tests/benchmark_kron_vs_direct.m')

clear; clc;

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'bvartools'));

if isOctave() == 0
    rng(123, 'twister');
else
    randn('state', 123);
end

% [ny, lags] pairs: small (~matches the test suite), medium (a typical
% empirical VAR), large (a sizeable system / small FAVAR-style model).
configs = [ ...
    3,  2; ...
    8,  4; ...
    15, 4  ...
    ];

n_reps = 200;   % repetitions per config, to get a stable timing average

fprintf('\n%-6s %-6s %-16s %-16s %-10s\n', 'ny', 'lags', 'old (ms/draw)', 'new (ms/draw)', 'speedup');
fprintf('%s\n', repmat('-', 1, 58));

for c = 1:size(configs, 1)
    ny   = configs(c, 1);
    lags = configs(c, 2);
    nk   = ny*lags + 1;

    % Valid lower-Cholesky factors of the right size. The actual entries
    % don't matter for timing, only the dimensions do, but they must be
    % genuine Cholesky factors (not arbitrary matrices) for the comparison
    % to reflect realistic floating point behavior.
    A1 = 0.1*randn(ny, ny);
    A2 = 0.1*randn(nk, nk);
    Sigma_lower_chol = chol(eye(ny) + A1*A1', 'lower');
    XXi_lower_chol   = chol(eye(nk) + A2*A2', 'lower');

    old_way = @() reshape(kron(Sigma_lower_chol, XXi_lower_chol) * randn(nk*ny, 1), nk, ny);
    new_way = @() XXi_lower_chol * reshape(randn(nk*ny, 1), nk, ny) * Sigma_lower_chol';

    t0 = tic;
    for r = 1:n_reps
        old_way();
    end
    t_old_ms = toc(t0) / n_reps * 1000;

    t0 = tic;
    for r = 1:n_reps
        new_way();
    end
    t_new_ms = toc(t0) / n_reps * 1000;

    fprintf('%-6d %-6d %-16.4f %-16.4f %-9.1fx\n', ny, lags, t_old_ms, t_new_ms, t_old_ms/t_new_ms);
end

fprintf('\nNote: this isolates only the Phi-draw step. A full bvar_ call also\n');
fprintf('spends time on OLS, IRFs, forecasts, etc., so the end-to-end speedup on\n');
fprintf('a real run will be smaller than the ratios above. The Phi-draw share of\n');
fprintf('total runtime -- and so the end-to-end benefit -- grows with ny, lags,\n');
fprintf('and the number of draws K.\n\n');
