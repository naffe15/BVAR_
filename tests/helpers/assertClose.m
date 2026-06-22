function [passed, message] = assertClose(actual, expected, tol, label)
%ASSERTCLOSE Compare two numeric arrays with a relative-tolerance check.
%
%   [passed, message] = assertClose(actual, expected, tol, label)
%
% Inputs:
%   actual, expected - numeric arrays of the same size
%   tol              - relative tolerance: elementwise
%                       |actual-expected| <= tol * max(1, |expected|)
%   label             - short text used in the returned message
%
% Outputs:
%   passed  - true/false
%   message - human-readable summary (used by run_tests.m instead of
%             throwing, so one failing case does not stop the rest of
%             the suite)
%
% Note on tolerance: a relative tolerance with a floor of 1 (rather than
% a pure relative test) avoids spurious failures on values that are
% close to zero, where a tiny absolute difference can look like a huge
% relative one.

if nargin < 4
    label = 'value';
end

if ~isequal(size(actual), size(expected))
    passed  = false;
    message = sprintf('%s: size mismatch (actual %s vs expected %s)', ...
        label, mat2str(size(actual)), mat2str(size(expected)));
    return
end

diffs  = abs(actual - expected);
scale  = max(1, abs(expected));
relerr = diffs ./ scale;
maxerr = max(relerr(:));

passed = maxerr <= tol;
if passed
    message = sprintf('%s: OK (max rel. error %.3g <= tol %.3g)', label, maxerr, tol);
else
    message = sprintf('%s: FAILED (max rel. error %.3g > tol %.3g)', label, maxerr, tol);
end

end
