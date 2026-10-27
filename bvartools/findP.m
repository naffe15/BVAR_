function P = findP(C,B,Q,p,k,index)
%FINDP Rotation that satisfies the short and long-run zero restrictions.
%
%   P = FINDP(C,B,Q,p,k,index)
%
%   Inputs
%     C       k x k     candidate impact matrix (from GENERATEDRAW)
%     B       (k*p+1) x k   VAR coefficients with the intercept in the
%                           first row, then the lag blocks
%     Q       k x 1 cell    restriction matrices from FINDQS
%     p       scalar    number of lags
%     k       scalar    number of variables
%     index   k x 1     original order of the shocks, from FINDQS
%
%   Output
%     P       k x k     orthogonal matrix; C*P satisfies the zero restrictions
%
%   See also FINDQS, GENERATEDRAW, IRESPONSE_ZEROS_SIGNS.
%
%   Reference: Rubio-Ramirez, J. F., D. F. Waggoner and T. Zha (2010),
%   Review of Economic Studies 77(2), 665-696.

%   Copyright Andrew Binning 2013. Please feel free to use and modify this
%   code as you see fit. If you use it in academic work, please cite
%   Binning, A. (2013), "Underidentified SVAR models: A framework for
%   combining short and long-run restrictions with sign-restrictions",
%   Norges Bank Working Paper 2013/14.

L0 = C;

beta_temp = B(2:end,:)';

beta = zeros(k,k);

for ii = 1:p
    
    beta = beta + beta_temp(:,(1:k)+(ii-1)*k);
    
end

Linf = (eye(k)-beta)\C;

F = [L0;Linf];

P = zeros(k,k);

for ii = 1:k
    if ii == 1
        Qtilde = Q{ii}*F;
    else
        Qtilde = [Q{ii}*F;P'];
    end
    [QQ,RR] = qr(Qtilde');
    P_temp = QQ(:,end);
    P(:,ii) = P_temp;
    
end

P = P(:,index);

