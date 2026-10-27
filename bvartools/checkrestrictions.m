function d = checkrestrictions(restriction,y,v)
%CHECKRESTRICTIONS Checks whether all restrictions in a list hold.
%
%   d = CHECKRESTRICTIONS(restriction,y)
%   d = CHECKRESTRICTIONS(restriction,y,v)
%
%   Inputs
%     restriction  cell array of strings, each a MATLAB expression in y
%                  and/or v that returns a logical scalar or array, e.g.
%                  'y(1,1:4,2)>0' or 'v(10,1)>0'
%     y            ny x hor x nshocks   impulse responses (can be [])
%     v            T x nshocks          structural shocks (default [])
%
%   Output
%     d            1 if every expression is true in all its elements, 0
%                  otherwise
%
%   See also CHECKRESTRICTIONS2, IRESPONSE_SIGN_NARRATIVE.

if nargin<3
    v = [];
end

d=0; 
count   = 0;
for ii = 1 : size(restriction,2)
    tmp = eval(restriction{ii});
    count = count + min(tmp);
end
if isempty(count)==1  
    error('There is a nan in the narrative restrictions.');
end
if count == size(restriction,2) % if all signs are verified stop
    d=1;
end
