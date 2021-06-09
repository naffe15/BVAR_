function J=jacob_bvar(param)

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018

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
