function Qbar  = max_fevd(i, h, j, Phi, Sigma, Kappa)
% finds the rotation where shock j maximizes the fevd of vriable i at horizon h  
% see also fevd.m

if nargin < 6
    Kappa = 1000;
end

N    = size(Sigma,1);
crit = nan(Kappa,1);
Q    = nan(N,N,Kappa);
for k = 1 : Kappa
    Q(:,:,k) = generateQ(N);                % generate an orthonormal matrix
    FEVD     = fevd(h,Phi,Sigma,Q(:,:,k));  % compute the FEVD
    % Calculate the contribution of shock j, to variable i forecast error volatility, at horizon h
    crit(k,1) = FEVD(i, h, j);
end
% Pick the maximum
[~,index] = max(crit);
Qbar      = Q(:,:,index); 

end
