function [ywc1, ywt1, ywc2, ywt2] = wavefilter(X, fig)

% wavelets (Haar) filtering: see Lubik, Matthes, Verona (2019)

%  ywc1 has the  cycles at frequencies 8-32 quarters (loose 16 datapoints)
%  ywc2 has the  cycles at frequencies 8-64 quarters (loose 32 datapoints)
%  ywt1 and  ywt2 have  the trends  and are defined  as  residuals
 

enddT=size(X,1);
ywc1=zeros(enddT,size(X,2));   ywt1=zeros(enddT,size(X,2)); 
ywc2=zeros(enddT,size(X,2));   ywt2=zeros(enddT,size(X,2));    

for qq=1:size(X,2)
%    y=zeros(enddT,1); 
    xx1=zeros(enddT,1); xx2=zeros(enddT,1); xx3=zeros(enddT,1);
    xx4=zeros(enddT,1); xx5=zeros(enddT,1); 

    y=squeeze(X(:,qq));

    % 2-4  quarters cycles
    for  tt=2:length(y)
    xx1(tt,1)=(1/2)*(y(tt)-(y(tt)-y(tt-1)));
    end
    % 4-8  quarters cycles
    for  tt=4:length(y)
    xx2(tt,1)=(1/4)*(y(tt)+y(tt-1)-(y(tt-2)+y(tt-3)) );
    end
    % 8-16  quarters cycles
    for  tt=8:length(y)
    xx3(tt,1)=(1/8)*(y(tt)+y(tt-1)+y(tt-2)+y(tt-3)- ...
                  (y(tt-4)+y(tt-5)+y(tt-6)+y(tt-7)) );
    end
    % 16-32 quarters cycles
    for  tt=16:length(y)
    xx4(tt,1)=(1/16)*(y(tt)+y(tt-1)+y(tt-2)+y(tt-3)+y(tt-4)+ ...
                    y(tt-5)+y(tt-6)+y(tt-7)-...
                   (y(tt-8)+y(tt-9)+y(tt-10)+y(tt-11)+y(tt-12)+ ...
                    y(tt-13)+y(tt-14)+y(tt-15)) );
    end
    %  32-64 quarters cycles
    for  tt=32:length(y)
     xx5(tt,1)=(1/32)*(y(tt)+y(tt-1)+y(tt-2)+y(tt-3)+y(tt-4)+y(tt-5)+ ...
        y(tt-6)+y(tt-7)+y(tt-8)+y(tt-9)+y(tt-10)+ y(tt-11)+y(tt-12)+ ...
         y(tt-13)+ y(tt-14)+y(tt-15) - ...
        (y(tt-16)+y(tt-17)+y(tt-18)+y(tt-19)+y(tt-20)+y(tt-21)+ ...
         y(tt-22)+y(tt-23)+y(tt-24)+y(tt-25)+y(tt-26)+y(tt-27)+ ...
         y(tt-28)+y(tt-29)+y(tt-30)+y(tt-31) ) );
    end
   % cyclical components 
   ywc1(16:length(y),qq)=xx3(16:length(y))+xx4(16:length(y));
   ywc2(32:length(y),qq)=xx3(32:length(y))+...
                       xx4(32:length(y))+xx5(32:length(y));
   % trend components (actual-BC cycles-high frequnecy cycles)
   ywt1(16:length(y),qq)=y(16:length(y))-ywc1(16:length(y),qq) ...
                         -xx1(16:length(y))-xx2(16:length(y));
   ywt2(32:length(y),qq)=y(32:length(y))-ywc2(32:length(y),qq) ...
                         -xx1(32:length(y))-xx2(32:length(y));
   
% xxx=[xx1 xx2 xx3 xx4  xx5]
% figure(100)
%  plot(xxx)
 
  if fig==1
      % plot  cyclical components
   figure(1)
   plot(ywc1(1:length(y),qq),'r', 'linewidth',2); hold  on; 
   plot(ywc2(1:length(y), qq),'b','linewidth',2); hold  off; axis  tight;
   legend('BC(8-32)','BC+LOW(8-64)')
   pause
   end
 end
 
 
end

