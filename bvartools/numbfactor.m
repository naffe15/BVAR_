function [A] = numbfactor(XX,TRANF)

%% Selecting the number of factor
%% Codes based on Bai and Ng (2002,ECMA)
%%

disp('                            ');
disp('Select the number of factors');

if nargin < 2
    prompt = 'Do you want to use raw (0), demeaned (1) or standardized (2) data ?';
    str2 = input(prompt,'s');
    TRANF   = str2num(str2);
end


kmax    = 20;%min(size(XX));
% colors  = strvcat('b','r','k','g','y');

for jj = 1 : 3
    [res] = baing(XX,kmax,jj,TRANF);
    A(jj) = res.ic1; 
%     plot(res.IC1,colors(jj,:));
%     hold on;
end
disp(sprintf('Demean %d',TRANF));
disp(sprintf('T = %d N = %d',size(XX)));
disp(sprintf('PCp1 = %d PCp2 = %d PCp3 = %d',A))

% disp('                            ');
% disp(sprintf('The min of all the Bai and Ng Criteria = %d', round(min(A))));
% nfac = round(min(A));
% disp(sprintf('Number of Factor = %d',nfac));

% if exist('nfac')==0,
%     prompt = 'How many factor you want to consider ?';
%     strf = input(prompt,'s');
%     nfac = str2num(strf);
% end

% if TRANF == 0, [ehat,fhat,lambda,ss]  = pc_T(XX,nfac);
% elseif TRANF == 1, [ehat,fhat,lambda,ss]  = pc_T(demean(XX),nfac);
% elseif TRANF == 2, [ehat,fhat,lambda,ss]  = pc_T(standard(XX),nfac);
% end
