function [nu] = u2nu(e_draws,Sigma_draws,Omega_draws)

[n,m,k] = size(Sigma_draws);

if nargin < 3
    Omega_draws = repmat(eye(n),1,1,k);
end

nu = zeros(size(e_draws));
for dd = 1 : k
    errors     = e_draws(:,:,dd);
    A          = chol( Sigma_draws(:,:,dd),'lower');
    nu(:,:,dd) = errors * inv( Omega_draws(:,:,dd)' * A'); %#ok<MINV>

    % nu(:,:,dd) = errors / (  Omega_draws(:,:,dd)' * A'); % structural innovations
end

