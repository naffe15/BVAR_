function [fis,fpt]= isdate(y,turnphase)
%  [fis,fpt]= isdate(y)
%
%  Determines if a data point is a potential turning point


adf=zeros(turnphase*2,1);


i=turnphase+1;


fis=0;
fpt=0;


q=1;
qq=1;
while q<=turnphase*2+1
    if q~=i        
        adf(qq) = indicat(y(q)-y(i));
        qq=qq+1;
    end
    q=q+1;
end
p1=sum(adf);

if p1==turnphase*2
    fis=1;
    fpt=1;
end

q=1; qq=1;
while q<=turnphase*2+1
    if q~=i        
        adf(qq) = indicat(y(i)-y(q));
        qq=qq+1;
    end
    q=q+1;
end
p1=sum(adf);
if p1==turnphase*2
    fis=1;
    fpt=-1;
end


function [ff] = indicat(ddy)
    ff=0;
    if ddy<=0
        ff=1;
    end
end

% End of function
end






