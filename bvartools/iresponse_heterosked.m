function [ir,Omeg] = iresponse_heterosked(Phi,errors,hor,delta,cont)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'iresponse_heterosked' computes the impulse response functions with
% Volatility-identified SVAR 
% Reference: Sims, C. A. (2020): "SVAR Identification Through Heteroskedasticity
% with Misspecified Regimes," Discussion paper, Princeton University,
% http://sims.princeton.edu/yftp/bpss/IDHmsspcfdRgms.pdf.

% Inputs:
% - errors, VAR reduced form errors (T x n)
% - Phi, AR parameters of the VAR
% - hor, horizon of the IRF
% - delta, binary variable with the regimes

% Output:
% - ir contains the IRF
% 1st dimension:   variable
% 2st dimension:   horizon
% 3st dimension:   shock

% Filippo Ferroni, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n]   = size(Phi);
favar   = 1;
if nargin < 5
    cont  = eye(n);
    favar = 0;
end
Omeg    = nan(n);
if favar == 0
    ir      = nan(n,hor,n);
else
    ir      = nan(size(cont,1),hor,n);
    yy      = nan(size(cont,1),hor,n);   
end
do_pagemtimes = 0;
if exist('pagemtimes','builtin') == 5
    do_pagemtimes = 1;
end

% compute the impact matrix 
S0      = cov(errors((delta==0),:));
S1      = cov(errors((delta==1),:));
X       = S0\S1;  %iS0S1 = inv(S0')*S1;
[~,S,V] = svd(X);
Omega   = V;
% compute impulse responses
y = iresponse(Phi,S,hor,Omega);
% uncompress the factors (if favar)
if favar == 1
    if do_pagemtimes ==0
        for jj  = 1 : size(y,3)
            yy(:,:,jj) = cont * y(:,:,jj);
        end
    else
        lam  = repmat(cont,1,1,n);
        yy   = pagemtimes(lam,y);
    end
    clear y; y = yy;
end
ir = y;
Omeg = Omega;