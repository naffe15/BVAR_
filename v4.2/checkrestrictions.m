function d = checkrestrictions(restriction,y,v)

% Check the restrictions 

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
