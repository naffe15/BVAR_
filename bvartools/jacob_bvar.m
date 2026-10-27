function J=jacob_bvar(param)
%JACOB_BVAR Jacobian of the exp transformation of the hyperparameters.
%
%   J = JACOB_BVAR(param)
%
%   BVAR_MAX_HYPER optimises over log hyperparameters; J converts the
%   Hessian back to their levels.
%
%   Input
%     param   n x 1   log hyperparameters
%
%   Output
%     J       n x n   diag(exp(param))
%
%   See also BVAR_MAX_HYPER.

%   Filippo Ferroni, 6/1/2015.

J=zeros(length(param),1);
for jj = 1 : length(param)
    J(jj)= bound0prime(param(jj));
end
J=diag(J);
 
% 
% function y = bound01prime(x);
% y = exp(x)/(1+exp(x))^2;
% 
function y = bound0prime(x);
y = exp(x);
