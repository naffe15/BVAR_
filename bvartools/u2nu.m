function [nu] = u2nu(e_draws,Sigma_draws,Omega_draws)

nu=zeros(size(e_draws));
for dd = 1 : size(e_draws,3)
    errors     = e_draws(:,:,dd);
    A          = chol( Sigma_draws(:,:,dd),'lower');
    nu(:,:,dd) = errors * inv( Omega_draws(:,:,dd)' * A'); %#ok<MINV>

    % nu(:,:,dd) = errors / (  Omega_draws(:,:,dd)' * A'); % structural innovations
end

