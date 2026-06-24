function run_tests()
%RUN_TESTS Master test runner for the BVAR_ toolbox.
%
% Usage:
%   run('<path to this file>/run_tests.m')
% or, with this folder on the path:
%   run_tests()
%
% Discovers every test_*.m function under cases/, runs each one, and
% prints a pass/fail summary table. Throws an error if any test failed
% (non-zero exit when run with -batch, so this can be wired into CI
% later without modification).
%
% Each test case is a plain function (no test-framework class
% dependency, deliberately -- this needs to run identically under both
% MATLAB and Octave) that takes no inputs and returns a struct with
% fields: name, passed, message.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'bvartools'));
addpath(fullfile(here, '..', 'cmintools'));
addpath(fullfile(here, 'helpers'));
addpath(fullfile(here, 'cases'));

case_dir = fullfile(here, 'cases');
files    = dir(fullfile(case_dir, 'test_*.m'));

if isempty(files)
    fprintf('No test cases found in %s\n', case_dir);
    return
end

n_pass = 0;
n_fail = 0;

fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('Running %d test case(s)\n', numel(files));
fprintf('%s\n\n', repmat('=', 1, 70));

for i = 1:numel(files)
    [~, fname] = fileparts(files(i).name);
    test_fun   = str2func(fname);
    try
        result = test_fun();
    catch err
        result.name    = fname;
        result.passed  = false;
        result.message = ['ERROR: ' err.message];
    end

    if result.passed
        status = 'PASS';
        n_pass = n_pass + 1;
    else
        status = 'FAIL';
        n_fail = n_fail + 1;
    end
    fprintf('[%s] %-40s %s\n', status, result.name, result.message);
end

fprintf('\n%s\n', repmat('-', 1, 70));
fprintf('%d passed, %d failed (of %d)\n', n_pass, n_fail, numel(files));
fprintf('%s\n\n', repmat('-', 1, 70));

if n_fail > 0
    error('run_tests:failures', '%d test case(s) failed.', n_fail);
end

end
