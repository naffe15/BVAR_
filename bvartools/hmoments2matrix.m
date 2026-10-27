function [f] = hmoments2matrix(hmoments,ny)
%HMOMENTS2MATRIX Marks which higher moments options.hmoments restricts.
%
%   f = HMOMENTS2MATRIX(hmoments,ny)
%
%   Inputs
%     hmoments  cell array of strings 'hm(a,b,c) > threshold' (see
%               CHECKMRESTRICTIONS)
%     ny        scalar   number of variables (shocks)
%
%   Output
%     f         ny x 4 x 2   1 where hm(a,b,c) is restricted, NaN elsewhere:
%                            a shock, b estimator (1 sample moments, 2-4
%                            robust estimators), c moment (1 skewness,
%                            2 kurtosis)
%
%   See also CHECKMRESTRICTIONS, IRESPONSE_SIGN_HMOMENTS.

f = nan(ny,4,2);
for ell = 1 : length(hmoments)
    tmpstr   = hmoments{ell};    
    startstr = findstr(tmpstr,'('); %#ok<*FSTR>
    endstr  = findstr(tmpstr,')');
    eval(['f' tmpstr(startstr:endstr) '=1;']) 
end

% % J = nan;
% for kk = 1 : 2
%     [row, col] = find(f(:,:,kk)==1);
%     J{kk} = [row, col];
% end

