function [f] = hmoments2matrix(hmoments,ny)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set of higher moments - contains the moments of 
% 1st dimension:   shocks
% 2nd dimension:   moment estimator (1= based on sample 3rd and 4th
% moments, >1 robust estimators
% 3rd dimension:   moment order (1 = skewnewss, 2 = kurtosis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

