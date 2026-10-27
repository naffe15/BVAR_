function [nu] = u2nu(e_draws,Sigma_draws,Omega_draws)
%U2NU Structural shocks from reduced-form residuals, draw by draw.
%
%   nu = U2NU(e_draws,Sigma_draws)
%   nu = U2NU(e_draws,Sigma_draws,Omega_draws)
%
%   Inputs
%     e_draws       T x n x K   residuals, e.g. BVAR.e_draws
%     Sigma_draws   n x n x K   residual covariances, e.g. BVAR.Sigma_draws
%     Omega_draws   n x n x K   rotations of chol(Sigma,'lower'), e.g.
%                               BVAR.Omegas (default eye: recursive shocks)
%
%   Output
%     nu            T x n x K   structural shocks
%
%   See also HISTDECOMP, BVAR_.

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

