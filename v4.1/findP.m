% Copyright Andrew Binning 2013
% Please feel free to use and modify this code as you see if fit. If you
% use this code in any academic work, please cite 
% Andrew Binning, 2013.
% "Underidentified SVAR models: A framework for combining short and long-run restrictions with sign-restrictions,"
% Working Paper 2013/14, Norges Bank.
function P = findP(C,B,Q,p,k,index)
%==========================================================================
% finds the rotation matrix that satisfies the short and long run
% restrictions.  Based on Juan F. Rubio-Ramirez & Daniel F. Waggoner & Tao Zha, 2010.
% "Structural Vector Autoregressions: Theory of Identification and
% Algorithms for Inference," Review of Economic Studies, Oxford University
% Press, vol. 77(2), pages 665-696.
% 
% inputs:
% C = initial short run impact matrix, usually from a cholesky
% decomposition of the forecast error variance
% B = Matrix of coefficients (including intercept estimates)
% Q = A cell containing the linear restrictions for each columnn
% p = number of lags
% k = number of dependent variables
% index = original column ordering in the matrix of restrictions
%
% outputs:
% P = orthogonal rotation matrix 
%==========================================================================

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

