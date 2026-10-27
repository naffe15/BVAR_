function [fis,fpt]= isdate(y,turnphase)
%ISDATE Checks whether the centre of a window is a local peak or trough.
%
%   [fis,fpt] = ISDATE(y,turnphase)
%
%   Inputs
%     y          (2*turnphase+1) x 1   window of the series
%     turnphase  scalar                half-width of the window
%
%   Outputs
%     fis   1 if the centre is a turning point, 0 otherwise
%     fpt   1 peak, -1 trough, 0 neither
%
%   See also RAWALL, DATE_.

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






