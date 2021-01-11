function [bcp5,bct5,nbp,nbt]=rawall(y,turnphase,nd,phase,cycle,thresh)
%  [bcp5,bct5,nbp,nbt]=rawall(y)
%
%  Makes sure all restrictions are imposed and chooses turning points


bcp5=[];
bct5=[];

adf =zeros(4,1);
i  = turnphase + 1;


while i<=nd-turnphase
    
    [fis,fpt] = isdate(y(i-turnphase:i+turnphase),turnphase);
    if fis==1&&(fpt==1)
        bcp5=[bcp5;i];
        
        break;
    elseif fis==1&&(fpt==-1)
        bct5=[bct5;i];
        
        break;
    end
    i=i+1;
end



if rows(bcp5)+rows(bct5)==0
    nbp=0;nbt=0;bcp5=0;bct5=0;
    
else
    
    
    if fpt==1
        j=bcp5;ltp=fpt;
        i=j+1;
        jj=0;
    else
        j=bct5;ltp=fpt;
        jj=bct5;
        i=j+1;
    end
    
    
    
    
    
    while i<=nd-turnphase
        [fis,fpt]=isdate(y(i-turnphase:i+turnphase),turnphase);
        if rows(bcp5)+rows(bct5)<2;
            if fis==1&(fpt==1)&(i-j>=phase)&(ltp==-1)&(indicat(y(i)-y(jj))==0);
                bcp5=[bcp5;i];
                ltp=fpt;
                j=i;
            elseif fis==1&(fpt==-1)&(i-j>=phase)&(ltp==1)&(indicat(y(j)-y(i))==0)|(fis==1)&(fpt==-1)&(y(i)-y(j)<-thresh)&(ltp==1);
                bct5=[bct5;i];
                ltp=fpt;
                jj=i;
                
            elseif fis==1&(fpt==1)&(ltp==1)&(indicat(y(j)-y(i))==1);
                
                bcp5 = [];
                bcp5=[bcp5;i];
                ltp=fpt;
                j=i;
                
            elseif fis==1&(fpt==-1)&(ltp==-1)&(indicat(y(i)-y(j))==1);
                bct5=[];
                bct5=[bct5;i];
                ltp=fpt;
                jj=i;
            end
            
        else;
            
            
            if fis==1&(fpt==1)&(i-jj>=phase)&(ltp==-1)&(i-j>=cycle)&(indicat(y(i)-y(jj))==0);
                bcp5=[bcp5;i];
                ltp=fpt;
                j=i;
                
                i=i+1;
                continue;
                
            elseif fis==1&(fpt==1)&(ltp==-1)&(indicat(y(i)-y(jj))==0);
                if rows(bct5)>1;
                    
                    if(y(i)<y(j))&(i-jj<phase)|(y(i)<y(j))&(i-jj<phase)&(i-j<cycle)|(y(i)<y(j))&(i-j<cycle)...
                            |((y(i)>=y(j))&(y(jj)<=y(bct5(rows(bct5)-1)))&(i-jj<phase)) ...
                            |((y(i)>=y(j))&(y(jj)<=y(bct5(rows(bct5)-1)))&(i-jj<phase)&(i-j<cycle));
                        
                    elseif(y(i)>=y(j))&(y(jj)>y(bct5(rows(bct5)-1)))&(i-jj<phase)|((y(i)>=y(j))&(y(jj)>y(bct5(rows(bct5)-1)))&(i-j<cycle));
                        bct5(rows(bct5),:)=[];
                        if rows(bcp5)==1;
                            bcp5=[];
                        else;
                            bcp5(rows(bcp5),:)=[];
                        end
                        bcp5=[bcp5;i];
                        ltp=fpt;
                        j=i;
                    elseif(y(i)>=y(j))&(y(jj)<=y(bct5(rows(bct5)-1)))&(i-j<cycle);
                        if rows(bct5)==2;
                            
                            bct5(1,:)=[];
                        else;
                            
                            
                            bct5t = bct5(rows(bct5));
                            bct5((rows(bct5)-1):rows(bct5),:) = [];
                            bct5 = [bct5;bct5t];
                            
                        end
                        if rows(bcp5)==1;
                            bcp5=[];
                        else;
                            bcp5(rows(bcp5),:)=[];
                        end
                        bcp5=[bcp5;i];
                        ltp=fpt;
                        j=i;
                    end
                elseif rows(bct5)==1;
                    if(y(i)<y(j))&(i-jj<phase)|(y(i)<y(j))&(i-jj<phase)&(i-j<cycle)|(y(i)<y(j))&(i-j<cycle)|((y(i)>=y(j))&(i-jj<phase)) ...
                            |((y(i)>=y(j))&(i-jj<phase)&(i-j<cycle));
                    elseif(y(i)>=y(j))&(i-j<cycle);
                        bcp5=[];
                        bcp5=[bcp5;i];
                        ltp=fpt;
                        j=i;
                    end
                    
                end
            elseif fis==1&(fpt==-1)&(i-j>=phase)&(ltp==1)&(i-jj>=cycle)&(indicat(y(j)-y(i))==0)|(fis==1)&(fpt==-1)&(y(i)-y(j)<-thresh)&(ltp==1);
                bct5=[bct5;i];
                ltp=fpt;
                jj=i;
                
                i=i+1;
                continue;
                
            elseif fis==1&(fpt==-1)&(ltp==1);
                
                if rows(bcp5) > 1;
                    if(y(i)>y(jj))&(i-j<phase)|(y(i)>y(jj))&(i-j<phase)&(i-jj<cycle)|(y(i)>y(jj))&(i-jj<cycle) ...
                            |((y(i)<=y(jj))&(y(j)>=y(bcp5(rows(bcp5)-1)))&(i-j<phase))       ...
                            |((y(i)<=y(j))&(y(j)>=y(bcp5(rows(bcp5)-1)))&(i-j<phase)&(i-jj<cycle));
                        
                    elseif(y(i)<=y(jj))&(y(j)<y(bcp5(rows(bcp5)-1)))&(i-j<phase)|((y(i)<=y(j))&(y(j)<y(bcp5(rows(bcp5)-1)))&(i-jj<cycle));
                        bcp5(rows(bcp5),:)=[];
                        if rows(bct5)==1;
                            bct5=[];
                        else;
                            bct5(rows(bct5),:)=[];
                        end
                        bct5=[bct5;i];
                        ltp=fpt;
                        jj=i;
                    elseif(y(i)<=y(jj))&(y(j)>=y(bcp5(rows(bcp5)-1)))&(i-jj<cycle);
                        if rows(bcp5)==2;
                            
                            bcp5(1,:)=[];
                        else;
                            
                            
                            bcp5t = bcp5(rows(bcp5));
                            bcp5((rows(bcp5)-1):rows(bcp5),:) = [];
                            bcp5 = [bcp5;bcp5t];
                            
                            
                        end
                        if rows(bct5)==1;
                            bct5=[];
                        else;
                            bct5(rows(bct5),:)=[];
                        end
                        bct5=[bct5;i];
                        ltp=fpt;
                        jj=i;
                    end
                    
                elseif rows(bcp5)==1;
                    
                    if(y(i)>y(jj))&(i-j<phase)|(y(i)>y(jj))&(i-j<phase)&(i-jj<cycle) ...
                            |(y(i)>y(jj))&(i-jj<cycle)...
                            |((y(i)<=y(jj))&(i-j<phase)) |(y(i)<=y(j))&(i-j<phase)&(i-jj<cycle);
                    elseif(y(i)<=y(jj))&(i-jj<cycle);
                        bct5=[];
                        bct5=[bct5;i];
                        ltp=fpt;
                        jj=i;
                        
                    end
                    
                end
                
            elseif(fis==1)&(fpt==1)&(ltp==1)&(indicat(y(j)-y(i))==1);
                if rows(bcp5)==1;
                    bcp5=[];
                else;
                    
                    bcp5(rows(bcp5),:)=[];
                end
                bcp5=[bcp5;i];
                ltp=fpt;
                j=i;
                
            elseif fis==1&(fpt==-1)&(ltp==-1)&(indicat(y(i)-y(jj))==1);
                if rows(bct5)==1;
                    bct5=[];
                else;
                    
                    bct5(rows(bct5),:)=[];
                end
                bct5=[bct5;i];
                ltp=fpt;
                jj=i;
            end
        end
        
        i=i+1;
    end
    
    
    
    
    if rows(bcp5)+rows(bct5)==1;
        nbp=0;nbt=0;bcp5=0;bct5=0;
    else;
        
        
        if (bcp5(1,:) > bct5(1,:));
            
            j = bct5(1,:)-1;
            while j > 1;
                if y(j,:)<y(bct5(1,:),:);
                    if rows(bct5)==1;
                        bct5 =[];
                    else;
                        bct5(1,:)=[];
                    end
                    break;
                end
                
                j = j-1;
            end
        elseif (bcp5(1,:) < bct5(1,:));
            
            j = bcp5(1,:)-1;
            while j > 1;
                if y(j,:)>y(bcp5(1,:),:);
                    if rows(bcp5)==1;
                        bcp5 =[];
                    else;
                        bcp5(1,:)=[];
                    end
                    
                    break;
                end
                
                j = j-1;
            end
            
        end
        
    end
    
    
    if rows(bcp5)+rows(bct5)==1;
        nbp=0;nbt=0;bcp5=0;bct5=0;
    else;
        
        
        nbp=rows(bcp5);nbt=rows(bct5);
        
        if (bcp5(nbp,:) > bct5(nbt,:));
            
            j = bcp5(nbp,:)+1;
            while j < nd;
                if y(j,:)>y(bcp5(nbp,:),:);
                    
                    if rows(bcp5)>2;
                        bcp5(rows(bcp5),:)=[];
                    else;
                        bcp5=[];
                    end
                    break;
                end
                
                j = j+1;
            end
            
        elseif (bcp5(nbp,:) < bct5(nbt,:));
            
            j = bct5(nbt,:)+1;
            while j < nd;
                if y(j,:)<y(bct5(nbt,:),:);
                    
                    if rows(bct5)>2;
                        bct5(rows(bct5),:)=[];
                    else;
                        bct5=[];
                    end
                    
                    break;
                end
                
                j = j+1;
            end
            
        end
        
    end
    
    nbp=rows(bcp5);nbt=rows(bct5);
    
end

function y = rows(x)
% Returns the row dimension of x
% USE:  y = rows(x)
[y,~] = size(x);
end

% function y = cols(x)
% % Returns the column dimension of x
% % USE:  y = cols(x)
% [~,y] = size(x);
% end
    
function [ff] = indicat(ddy)
    ff=0;
    if ddy<=0
        ff=1;
    end
end

% End of function
end