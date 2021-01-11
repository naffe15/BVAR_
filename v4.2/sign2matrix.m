function [f,sr] = sign2matrix(signs,ny)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yr  = ones(ny);
ys  = ones(ny);
y   = nan(ny);
for ell = 1 : length(signs)
    eval(signs{ell})    
end

f(1:ny , 1:ny)       = ys;
f(1+ny:2*ny , 1:ny)  = yr;
sr                   = y;


