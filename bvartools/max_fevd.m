function Qbar  = max_fevd(i, h, j, Phi, Sigma, Kappa)
% finds the rotation where shock j maximizes the fevd of vriable i at horizon h  
% see also fevd.m

if nargin < 6
    Kappa = 1000;
end
if length(i) ~= length(j)
    error('each shock must have one target')
end

N    = size(Sigma,1);
crit = nan(Kappa,length(j));
Q    = nan(N,N,Kappa);
for k = 1 : Kappa
    Q(:,:,k) = generateQ(N);                % generate an orthonormal matrix
    FEVD     = fevd(h,Phi,Sigma,Q(:,:,k));  % compute the FEVD
    % Calculate the contribution of shock j, to variable i forecast error volatility, at horizon h
    for j1 = 1 : length(j)
        crit(k,j1) = FEVD(i(j1), j(j1));
    end
end
% Pick the maximum
Qbar = nan(N,N,length(j));
for j1 = 1 : length(j)
    [~,index] = max(crit(:,j1));
    tmp_notindx = setdiff([1:N] , j(j1));
    Omega = Q(:,:,index);
    Qbar(:,:,j1) = Omega(:,[j(j1) tmp_notindx]);
end

end
