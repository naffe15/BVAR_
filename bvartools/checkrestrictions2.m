function [d,fsign] = checkrestrictions2(restriction,y)
%CHECKRESTRICTIONS2 Checks sign restrictions, allowing a sign flip of the shocks.
%
%   [d,fsign] = CHECKRESTRICTIONS2(restriction,y)
%
%   Inputs
%     restriction  cell array of strings, each a MATLAB expression in y,
%                  e.g. 'y(1,1:4,2)>0'
%     y            ny x hor x nshocks   impulse responses
%
%   Outputs
%     d            1 if all restrictions hold for y or for -y, 0 otherwise
%     fsign        1, or -1 when they hold for -y (the rotation must then
%                  be multiplied by -1)
%
%   See also CHECKRESTRICTIONS, IRESPONSE_SIGN.

d     = 0; 
count = 0;
fsign = 1;

for ii = 1 : size(restriction,2)
    tmp = eval(restriction{ii});
    count = count + min(tmp);
end
if isempty(count)==1  
    error('There is a nan in the narrative restrictions.');
end
% check if -1*IRF are verified
y        = -y;
count1   = 0;
for ii = 1 : size(restriction,2)
    tmp = eval(restriction{ii});
    count1 = count1 + min(tmp);
end
if count == size(restriction,2) || count1 == size(restriction,2) % if all signs are verified stop    
    d=1;
end
if count1 == size(restriction,2)
    fsign = -1;
end
