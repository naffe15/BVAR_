function FE = fe_(hor,errors,Phi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the in-sample forecast error at horizon h
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T , N]     = size(errors);
[m , k]     = size(Phi);
if rem(m, k)==0
    lags = m/k;
else
    lags = floor((m-1)/k);
end
MA = var2ma(Phi,hor); 

FE              = nan(T,N,hor);
FE(:,:,1)       = errors;
% FE(2:end,:,2)   = error(1:end-1,:)*MA(:,:,1)' +  error(2:end,:);
% FE(3:end,:,3)   = error(1:end-2,:)*MA(:,:,2)' +  error(2:end-1,:)*MA(:,:,1)' +  error(3:end,:);

for hh  = 2 : hor
    %if hh > lags, keyboard; end
    tmp_= errors(hh:end,:);
    for kk = 1 : hh-1
        tmp_ = tmp_ + errors(kk:end-(hh-kk),:)*MA(:,:,hh-kk)' ;
    end
    FE(hh:end,:,hh) = tmp_;
    clear tmp_
end



% % companion
% F       = [Phi(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
% G       = eye(N * lags, N);
% %C       = [Phi(end, :)'; zeros(N*(lags-1),1)];
% 
% A     = chol(Sigma,'lower');
% Kappa = G * A * Omega ;
% 
% tmp_=0;
% for hh  = 1 : hor
%     tmp_ = tmp_ + F^(hh-1) * Kappa * Kappa' * F^(hh-1)';
% end
% out_.all_var_  =  diag(tmp_(1:N,1:N));
% 
% for sho =1 : N
%     tmp_1 = 0;
%     Ind              = zeros(N);
%     Ind(sho,sho) = 1;
%     for hh  = 1 : hor
%         tmp_1 = tmp_1 + F^(hh-1) * Kappa * Ind * Kappa' * F^(hh-1)';
% 
%     end
%     out_.var_(:,sho) = diag(tmp_1(1:N,1:N));
% end
% 
% if max(max(abs( sum(out_.var_,2) - out_.all_var_))) > 1e-10
%     error('Something went wrong')
% end
% 
% for indx_sho = 1 : N
%     FEVD(:,indx_sho) = out_.var_(:,indx_sho)./out_.all_var_*100;    
% 
% end



