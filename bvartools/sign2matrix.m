function [f,sr] = sign2matrix(signs,ny)
%SIGN2MATRIX Converts options.zeros_signs strings into restriction matrices.
%
%   [f,sr] = SIGN2MATRIX(signs,ny)
%
%   Inputs
%     signs   cell array of strings, each evaluated as a MATLAB statement:
%               'y(j,k)=1;'  or 'y(j,k)=-1;'  sign of the impact of shock k
%                                             on variable j
%               'ys(j,k)=0;'                  zero impact (short run)
%               'yr(j,k)=0;'                  zero long-run effect
%     ny      scalar   number of variables
%
%   Outputs
%     f       2ny x ny   [ys; yr]: 0 marks a zero restriction, 1 no restriction
%     sr      ny x ny    impact signs: +1, -1, or NaN when unrestricted
%
%   See also IRESPONSE_ZEROS_SIGNS, FINDQS.

%   Filippo Ferroni, 6/1/2015.

yr  = ones(ny);
ys  = ones(ny);
y   = nan(ny);
for ell = 1 : length(signs)
    eval(signs{ell})    
end

f(1:ny , 1:ny)       = ys;
f(1+ny:2*ny , 1:ny)  = yr;
sr                   = y;


