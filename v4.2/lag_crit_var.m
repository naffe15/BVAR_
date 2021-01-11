function lag_crit_var(y,maxlag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'lag_crit_var' prints various info crit (not allowed NaN)

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

crt_ = nan(4,maxlag);

opt.K =1;
for ilag = 1 : maxlag
    tmp_ = bvar(y,ilag,opt);
    crt_(1,ilag) =  tmp_.InfoCrit.AIC;
    crt_(2,ilag) =  tmp_.InfoCrit.SIC;
    crt_(3,ilag) =  tmp_.InfoCrit.HQIC;
    crt_(4,ilag) =  tmp_.InfoCrit.BIC;
end
rownam = {'AIC','SIC','HQIC','BIC'};
for ilag = 1 : maxlag
    disp('=============================================')
    X = sprintf('%s = %0.0g','Number of lags',ilag);
    disp(X)
    for jj =1: size(crt_,1)
        X = sprintf('%s = %0.5g',rownam{jj},crt_(jj,ilag));
        disp(X)
    end
end
