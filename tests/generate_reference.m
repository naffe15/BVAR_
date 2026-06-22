%GENERATE_REFERENCE One-time script to (re)create the frozen "golden
% master" fingerprint used by tests/cases/test_bvar_minnesota_regression.m.
%
% Run this manually -- it is NOT part of run_tests.m and is not run
% automatically. Each run overwrites reference/bvar_minnesota_fingerprint.mat,
% which is the baseline that the regression test checks future runs
% against. Re-generate it only when you have manually verified that a
% change to bvar_ is an intentional, correct change -- never just to make
% a failing test pass, since that would defeat the point of the test.
%
% Usage (from any working directory):
%   run('<path to this file>/generate_reference.m')

close all; clear; clc;

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'bvartools'));
addpath(fullfile(here, '..', 'cmintools'));
addpath(fullfile(here, 'helpers'));

ny   = 3;
lags = 2;
T    = 300;
% Seed 7 here is deliberately different from the seed used in
% test_bvar_minnesota_recovery.m (42) -- this dataset only needs to be
% fixed and reproducible, not connected to a "known truth" check.
[y, ~, ~] = simulate_var(ny, lags, T, 7);

options              = struct();
options.priors.name  = 'Minnesota';
options.K            = 50;   % kept below 100 so bvar_ does not pop up a GUI waitbar
options.hor          = 4;
options.fhor         = 1;
options.noprint      = 1;

BVAR = bvar_(y, lags, options);

ir_mean = mean(BVAR.ir_draws, 4);   % dims: variable x horizon x shock x draw

fingerprint.Phi_mean   = mean(BVAR.Phi_draws, 3);
fingerprint.Sigma_mean = mean(BVAR.Sigma_draws, 3);
fingerprint.ir_h1      = squeeze(ir_mean(:,1,:));   % response at horizon 1, all var x shock pairs
fingerprint.logmlike   = BVAR.logmlike;

ref_dir = fullfile(here, 'reference');
if ~exist(ref_dir, 'dir')
    mkdir(ref_dir);
end
save(fullfile(ref_dir, 'bvar_minnesota_fingerprint.mat'), 'fingerprint');

fprintf('Reference fingerprint written to %s\n', ...
    fullfile(ref_dir, 'bvar_minnesota_fingerprint.mat'));
