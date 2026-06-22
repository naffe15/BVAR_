# BVAR_ test harness

A lightweight, MATLAB/Octave-compatible test suite for the toolbox. No
test-framework dependency (no `matlab.unittest` classes) so the same
code runs under both MATLAB and Octave.

## Layout

```
tests/
  run_tests.m              master runner: discovers and runs cases/test_*.m
  generate_reference.m     one-time script to (re)create golden-master fingerprints
  helpers/
    assertClose.m          relative-tolerance comparison used by all test cases
    simulate_var.m         generates synthetic VAR(p) data with a known, true Phi/Sigma
  cases/
    test_bvar_minnesota_recovery.m     correctness check (no stored reference)
    test_bvar_minnesota_regression.m   regression check (needs the reference file below)
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
