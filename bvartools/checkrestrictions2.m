function [d,fsign] = checkrestrictions2(restriction,y)

% Check the restrictions 
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
