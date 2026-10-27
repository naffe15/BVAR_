function [ir,Q]=iresponse_longrun(alpha,Sigma,hor,lags)
%IRESPONSE_LONGRUN Impulse responses with long-run (Blanchard-Quah) restrictions.
%
%   [ir,Q] = IRESPONSE_LONGRUN(alpha,Sigma,hor,lags)
%
%   The long-run impact matrix is lower triangular, so only the first shock
%   has a permanent effect on the first variable. The impact of the first
%   shock on the first variable is normalised to be positive.
%
%   Inputs
%     alpha   (ny*lags+nx) x ny  VAR coefficients, lag-1 block first
%     Sigma   ny x ny            covariance of the reduced-form residuals
%     hor     scalar             number of horizons (horizon 1 is the impact)
%     lags    scalar             number of lags
%
%   Outputs
%     ir      ny x hor x ny      responses (variable x horizon x shock)
%     Q       ny x ny            impact matrix
%
%   See also IRESPONSE, IRESPONSE_SET.
%
%   Reference: Blanchard, O. J. and D. Quah (1989), "The Dynamic Effects
%   of Aggregate Demand and Supply Disturbances", American Economic Review
%   79(4), 655-673.

%   Filippo Ferroni, 6/1/2015.

N            = size(Sigma,1);
ir          = zeros(N,hor,N); % variables, horizon, shock

% companion
F       = [alpha(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
G       = eye(N * lags, N);
Inp     = eye(size(F));
% long run sum of VMA 
C1             = Inp - F;
C1             = C1 \ Inp; % = inv(C1)*Inp;
C1             = C1(1 : N, 1 : N);
PSI1           = chol(C1 * Sigma * C1')';
% IMPACT MATRIX
Q              = (C1 \ eye(size(C1))) * PSI1;
% normalize to positive entry in the element Q(1,1)
Q(1,1)     = Q(1,1) * sign(Q(1,1));
Fk         = eye(N * lags);

for k=1:hor
    PHI         = Fk * G * Q;
    ir(1:N,k,:) = G'* PHI;
    Fk          = F * Fk;
end

    