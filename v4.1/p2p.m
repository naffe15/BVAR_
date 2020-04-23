function [posterior,prior] = p2p(hh, shrinkage,prior0,olsreg,F,G,Fo,positions_nylags,position_constant)
% 
%********************************************************
% Conjugate Prior: MN-IW for Direct Methods
%********************************************************

% settings
if isempty(position_constant) == 0
    nx = 1;
else 
    nx = 0;
end

ny       = size(G,2);
nylags   = size(F,1);
iIminusF = inv(eye(nylags) - F); 

% constructing the prior mean
Fhh        = F^hh;
Fohh       = (iIminusF * (eye(nylags)-Fhh)) * Fo;
if hh < 40
    priorSigmaScale  = zeros(ny);
    for j = 0 : hh
        priorSigmaScale = priorSigmaScale + G' * F^(hh-j)* G * prior0.Sigma.scale * G' * F^(hh-j)' * G;
    end
else
    % solves x-a*x*a'=b for b (and then x) symmetrical
    % function [x,info]=lyapunov_symm(a,b)
    [priorSigmaScale0,~] = lyapunov_symm(F,G*prior0.Sigma.scale*G' - F^(hh+1)* G * prior0.Sigma.scale * G' * F^(hh+1)');
    % max(max(abs(priorSigmaScale-G'*priorSigmaScale0*G)))
    priorSigmaScale      = G'*priorSigmaScale0*G;
end

prior.BetaMean   = [Fhh(1:ny,:)'; Fohh(1 : ny,1)'];
prior.BetaVar    = prior0.Phi.cov * 1/ shrinkage;

prior.df  = prior0.Sigma.df;            % usually number of regressors minus the 2
prior.XXi = inv( prior.BetaVar );       % V^{-1}
prior.S   = priorSigmaScale;            % Sigma0

% retrieve the OLS
B_   = olsreg.beta([positions_nylags position_constant], :);
XX_  = olsreg.X(:,[positions_nylags position_constant])' * olsreg.X(:,[positions_nylags position_constant]);
E_   = olsreg.error;
XXp_ = XX_ + prior.XXi;

% construct the posterior
posterior.df     = olsreg.N - nylags - nx + prior0.Sigma.df;
posterior.XXi    = inv( XXp_ );
posterior.PhiHat = posterior.XXi * (XX_ * B_ + prior.XXi * prior.BetaMean);
posterior.S     =  ...
    E_'* E_ + priorSigmaScale + prior.BetaMean' * prior.XXi * prior.BetaMean + ...
    B_' * XX_ * B_ - posterior.PhiHat' * XXp_ * posterior.PhiHat;
% FF    = prior.Sigma.scale + (posterior1.PhiHat - prior1.BetaMean)'* prior1.XXi * (posterior1.PhiHat - prior1.BetaMean) ...
%     + posterior1.E_'* posterior1.E_; 

posterior.B_   = B_ ;
posterior.XX_  = XX_;       %olsreg.X(:,[positions_nylags position_constant])' * olsreg.X(:,[positions_nylags position_constant]);
posterior.E_   = E_;        %olsreg.error';
posterior.XXp_ = XXp_;      %XX_ + prior.XXi;
posterior.U_   = olsreg.Y - olsreg.X(:,[positions_nylags position_constant]) * posterior.PhiHat;% posterior errors;
