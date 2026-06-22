# BVAR_ test harness

A lightweight, MATLAB/Octave-compatible test suite for the toolbox. No
test-framework dependency (no `matlab.unittest` classes) so the same
code runs under both MATLAB and Octave.

## Layout

```
tests/
  run_tests.m                  master runner: discovers and runs cases/test_*.m
  generate_reference.m         one-time script to (re)create golden-master fingerprints
  benchmark_kron_vs_direct.m   standalone, manual: timing comparison for the Phi-draw optimization
  helpers/
    assertClose.m          relative-tolerance comparison used by all test cases
    simulate_var.m         generates synthetic VAR(p) data with a known, true Phi/Sigma
  cases/
    test_bvar_minnesota_recovery.m     correctness check (no stored reference)
    test_bvar_minnesota_regression.m   regression check (needs the reference file below)
    test_bvar_minnesota_timing.m       coarse smoke test, generous threshold (see "Performance" below)
    test_bvar_penalized.m              smoke test for the Ridge/Lasso/ElasticNet path
  reference/
    bvar_minnesota_fingerprint.mat     created by generate_reference.m (not checked in yet)
```

## Two kinds of test, on purpose

**Recovery / correctness tests** (`*_recovery.m`) simulate data from a VAR
with *known* coefficients and check that the estimator's posterior mean
lands close to the truth. Loose tolerance since very low number of draws 
(50 iterations).  

**Regression tests** (`*_regression.m`) run a function on a fixed dataset
and compare the result to a frozen "golden master" saved earlier. They
catch *changes* in behavior, intended or not, but can't tell you whether
the original behavior was correct in the first place — that's what the
recovery tests are for.

`bvar_` resets its own RNG seed internally (`rng(999)` near the top of
`bvar_.m`), so repeated calls with identical inputs should be essentially
bit-for-bit reproducible. That's why the regression test uses a tight
tolerance (1e-6) while the recovery test uses a loose one (0.35) to
allow for prior shrinkage and small-sample noise.

There's also a third kind, **smoke tests** (`test_bvar_penalized.m`):
they just check that a code path runs and produces the right shape of
sane, finite output. No tolerance tuning involved — they exist to catch
"this errors now" or "this field is missing now", not "this number
changed by more than expected."

## Performance

`test_bvar_minnesota_timing.m` is a coarse smoke test on runtime, not a
real benchmark — its pass/fail threshold (60s) is deliberately loose
since wall-clock time depends on the machine. It exists to catch a gross
regression (10x+), not to certify a speed. It always reports the actual
measured time in its message, so it's worth glancing at over time even
when it passes.

For an actual measurement of the kron-vs-direct optimization in bvar_'s
Phi draw (see the comment in `bvar_.m` above the `Phi3 = XXi_lower_chol *
Z * Sigma_lower_chol'` line), run `benchmark_kron_vs_direct.m` manually.
It isolates just that one step and reports a speedup ratio across a few
problem sizes, which is far more informative than the timing test's
pass/fail outcome and doesn't depend on a chosen threshold at all.

## How to run

```matlab
run('/path/to/BVAR_/tests/run_tests.m')
```

The first time, `test_bvar_minnesota_regression` will fail because
`reference/bvar_minnesota_fingerprint.mat` doesn't exist yet. Create it
once with:

```matlab
run('/path/to/BVAR_/tests/generate_reference.m')
```

Then re-run `run_tests.m`. From then on, commit
`reference/bvar_minnesota_fingerprint.mat` to git so the regression test
is checking against a baseline everyone shares — re-generate it
deliberately (and only after manually verifying a change is correct),
never just to make a failing test pass.


## Extending

To add a test for another function, copy the pattern: write a
`test_<name>.m` in `cases/` that returns `struct('name', ..., 'passed',
..., 'message', ...)`, drop it in `cases/`, and `run_tests.m` will pick
it up automatically. Good next candidates, in rough priority order:
`cvar_` (deterministic OLS — exact-equality test, no RNG involved),
`iresponse_sign` (the most commonly used identification scheme),
`forecasts`/`cforecasts`, and `connectedness`.

`test_bvar_penalized.m` only smoke-tests Ridge unconditionally; it
skips Lasso/ElasticNet if the Statistics and Machine Learning Toolbox's
`lasso()` isn't on the path. If you regularly use those, worth checking
they actually ran (vs. silently skipped) by reading the printed message.
