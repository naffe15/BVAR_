function C = generateDraw(C,k)
%GENERATEDRAW Random impact matrix consistent with the residual covariance.
%
%   C = GENERATEDRAW(C,k)
%
%   Inputs
%     C   k x k    initial impact matrix, usually chol(Sigma,'lower')
%     k   scalar   number of variables
%
%   Output
%     C   k x k    C times a random orthonormal matrix (Haar distributed)
%
%   See also FINDP, FINDQS, GENERATEQ, IRESPONSE_ZEROS_SIGNS.
%
%   Reference: Rubio-Ramirez, J. F., D. F. Waggoner and T. Zha (2010),
%   Review of Economic Studies 77(2), 665-696.

%   Copyright Andrew Binning 2013. Please feel free to use and modify this
%   code as you see fit. If you use it in academic work, please cite
%   Binning, A. (2013), "Underidentified SVAR models: A framework for
%   combining short and long-run restrictions with sign-restrictions",
%   Norges Bank Working Paper 2013/14.

newmatrix = randn(k,k);

[Q,R] = qr(newmatrix);

for ii = 1:k
    if R(ii,ii)<0
        Q(:,ii) = -Q(:,ii);
    end
end

C = C*Q;