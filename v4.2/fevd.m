function FEVD = fevd(hor,Phi,Sigma,Omega)

% computes the forecast error variance decomposition

N           = size(Sigma,1);
[m , k]     = size(Phi);
lags        = (m-1)/k;

if nargin < 4
    Omega = eye(N);
end 


% companion
F       = [Phi(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
G       = eye(N * lags, N);
%C       = [Phi(end, :)'; zeros(N*(lags-1),1)];

A     = chol(Sigma,'lower');
Kappa = G * A * Omega ;

tmp_=0;
for hh  = 1 : hor
    tmp_ = tmp_ + F^(hh-1) * Kappa * Kappa' * F^(hh-1)';
end
out_.all_var_  =  diag(tmp_(1:N,1:N));

for sho =1 : N
    tmp_1 = 0;
    Ind              = zeros(N);
    Ind(sho,sho) = 1;
    for hh  = 1 : hor
        tmp_1 = tmp_1 + F^(hh-1) * Kappa * Ind * Kappa' * F^(hh-1)';
        
    end
    out_.var_(:,sho) = diag(tmp_1(1:N,1:N));
end

if max(max(abs( sum(out_.var_,2) - out_.all_var_))) > 1e-10
    error('Something went wrong')
end

for indx_sho = 1 : N
    FEVD(:,indx_sho) = out_.var_(:,indx_sho)./out_.all_var_*100;    

end



