function w = matrictint(S, df, XXi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the log of the integral of the kernel of the PDF of a
% normal-inverse-Wishart distribution.
%
% S:   parameter of inverse-Wishart distribution
% df:  number of degrees of freedom of inverse-Wishart distribution
% XXi: first component of VCV matrix of matrix-normal distribution
% 
% Computes the integral over (Phi, Sigma) of:
%
% det(Sigma)^(-k/2)*exp(-0.5*Tr((Phi-PhiHat)'*(XXi)^(-1)*(Phi-PhiHat)*Sigma^(-1)))*
% det(Sigma)^((df+ny+1)/2)*exp(-0.5*Tr(Sigma^(-1)*S))
%
% (where k is the dimension of XXi and ny is the dimension of S and
% Sigma)

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/matrictint.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
