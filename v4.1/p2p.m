function [posterior,prior] = p2p(hh,prior,olsreg,F,G,Fo,positions_nylags,position_constant)

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
iIminusF = inv(eye(nylags) - F); % iIminusF = (I-F)\eye(nylags);

% constructing the prior mean
Fhh        = F^hh;
Fohh       = (iIminusF * (eye(nylags)-Fhh)) * Fo;
if hh < 40
    priorSigmaScale  = zeros(ny);
    for j = 0 : hh
        priorSigmaScale = priorSigmaScale + G' * F^(hh-j)* G * prior.Sigma.scale * G' * F^(hh-j)' * G;
    end
else
    % solves x-a*x*a'=b for b (and then x) symmetrical
    % function [x,info]=lyapunov_symm(a,b)
    [priorSigmaScale0,~] = lyapunov_symm(F,G*prior.Sigma.scale*G' - F^(hh+1)* G * prior.Sigma.scale * G' * F^(hh+1)');
    % max(max(abs(priorSigmaScale-G'*priorSigmaScale0*G)))
    priorSigmaScale      = G'*priorSigmaScale0*G;
end

priorBetaMean   = [Fhh(1:ny,:)'; Fohh(1 : ny,1)'];
priorBetaVar    = prior.Phi.cov * 1/prior.tau(hh);
% Ghh        = G' * (iIminusF * (eye(nylags)-Fhh)) * G;
% priorSigmaScale = Ghh * prior.Sigma.scale * Ghh';

prior.df  = prior.Sigma.df; % usually number of regressors minus the 2
% prior.XXi = inv(prior.Phi.cov(hh * ny + 1 : hh * ny + 1 : ));
prior.XXi = inv( priorBetaVar );% \ one(nylags-1); % nendolags * nendo
prior.S   = priorSigmaScale;
% prior_int = matrictint(priorSigmaScale, prior.df, prior.XXi); % ny * ny

% retrieve the OLS
B_   = olsreg.beta([positions_nylags position_constant], :);
XX_  = olsreg.X(:,[positions_nylags position_constant])' * olsreg.X(:,[positions_nylags position_constant]);
E_   = olsreg.error';
XXp_ = XX_ + prior.XXi;

% construct the posterior
posterior.df     = olsreg.N - nylags - nx + prior.Sigma.df;
posterior.XXi    = inv( XXp_ );
posterior.PhiHat = posterior.XXi * (XX_ * B_ + prior.XXi * priorBetaMean);
posterior.S     =  ...
    E_'* E_ + priorSigmaScale + priorBetaMean' * prior.XXi * priorBetaMean + ...
    B_' * XX_ * B_ - posterior.PhiHat' * XXp_ * posterior.PhiHat;
