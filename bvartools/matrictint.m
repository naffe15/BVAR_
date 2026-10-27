function w = matrictint(S, df, XXi)
%MATRICTINT Log integral of a normal-inverse-Wishart kernel.
%
%   w = MATRICTINT(S,df,XXi)
%
%   Computes the log of the integral over (Phi,Sigma) of
%     det(Sigma)^(-k/2) * exp(-0.5*tr((Phi-PhiHat)'*inv(XXi)*(Phi-PhiHat)*inv(Sigma)))
%     * det(Sigma)^(-(df+ny+1)/2) * exp(-0.5*tr(inv(Sigma)*S))
%   where k is the size of XXi and ny the size of S.
%
%   Inputs
%     S     ny x ny   scale of the inverse-Wishart
%     df    scalar    degrees of freedom of the inverse-Wishart
%     XXi   k x k     first component of the covariance of the
%                     matrix-normal part
%
%   Output
%     w     scalar    log of the integral
%
%   See also BVAR_ML, MNIW_LOG_DNSTY.

%   Original file by Christopher A. Sims, downloaded from
%   http://sims.princeton.edu/yftp/VARtools/matlab/matrictint.m

k=size(XXi,1);
ny=size(S,1);
[cx,p]=chol(XXi);
[cs,q]=chol(S);

if any(diag(cx)<100*eps)
    error('singular XXi')
end
if any(diag(cs<100*eps))
    error('singular S')
end

% Matrix-normal component
w1 = 0.5*k*ny*log(2*pi)+ny*sum(log(diag(cx)));

% Inverse-Wishart component
w2 = -df*sum(log(diag(cs))) + 0.5*df*ny*log(2) + ny*(ny-1)*0.25*log(pi) + ggammaln(ny, df);

w = w1 + w2;

function lgg = ggammaln(m, df)
if df <= (m-1)
    error('Too few df in ggammaln: increase the # of obs or decrease the # of lags.')
else
    garg = 0.5*(df+(0:-1:1-m));
    lgg = sum(gammaln(garg));
end
