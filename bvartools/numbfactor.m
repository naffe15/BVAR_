function [A] = numbfactor(XX,TRANF)
%NUMBFACTOR Number of static factors by the Bai and Ng (2002) criteria.
%
%   A = NUMBFACTOR(XX)
%   A = NUMBFACTOR(XX,TRANF)
%
%   Inputs
%     XX      T x N      panel, time in rows
%     TRANF   0, 1 or 2  data treatment: 0 raw, 1 demeaned, 2 standardized.
%                        If omitted, the user is asked for it.
%
%   Output
%     A       1 x 3      number of factors selected by PCp1, PCp2 and PCp3,
%                        searching up to kmax = 20; also printed to the
%                        command window
%
%   Local functions: BAING (criterion values for one penalty), MININDC
%   (row index of the minimum of each column), SEQA (additive sequence).
%
%   See also PC_T, DEMEAN, STANDARD.
%
%   Reference: Bai, J. and S. Ng (2002), "Determining the Number of Factors
%   in Approximate Factor Models", Econometrica 70(1), 191-221.
%   Adapted from the code of Bai and Ng.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res]=baing(x,kmax,jj,DEMEAN)

% BAING computes the Bai-Ng criterion jj (1,2,3 = PCp1,PCp2,PCp3) for
% 0..kmax factors and returns the minimiser in res.ic1.
% X is observed
% r is the true number of true factors
% F is T by r matrix of true factors
% Lambda N by r is the true loading matrix
% C=F*Lambda' T by N is the true common component
% chat is the estimated common component
% DEMEAN=2 standardize data, DEMEAN=1, demean data, DEMEAN=0 raw data
  
T   = size(x,1);
N   = size(x,2);
NT  = N*T;
NT1 = N+T;
CT  = zeros(1,kmax);
ii  = 1 : 1 : kmax;

if jj ==1;
    CT(1,:) = log(NT/NT1)*ii*NT1/NT;    
elseif jj==2; 
    CT(1,:) = (NT1/NT)*2*log(min([N;T]))*ii;
elseif jj==3; 
    GCT     = min([N;T]);
    CT(1,:) = ii*log(GCT)/GCT; 
end;

if DEMEAN ==2;
    X=standard(x);
elseif DEMEAN ==1;
    X=demean(x);
elseif DEMEAN==0;
    X=x;
end;

IC1                  = zeros(size(CT,1),kmax+1);
Sigma                = zeros(1,kmax+1);
XX                   = X*X';
[Fhat0,eigval,Fhat1] = svd(XX);
for i = kmax : -1 : 1,
    % Fhat   = Fhat0(:,1:i);
    % lambda = Fhat'*X;
     
    fhat      = Fhat0(:,1:i)*sqrt(T);  
    lambda    = X'*fhat/T;
    chat      = fhat * lambda';
    ehat      = X - chat;
    Sigma(i)  = mean(sum(ehat.*ehat/T));
    IC1(:,i)  = log(Sigma(i)) + CT(:,i);
end
Sigma(kmax+1)  = mean(sum(X.*X/T));
IC1(:,kmax+1)  = log(Sigma(kmax+1));
ic1            = minindc(IC1')';
ic1            = ic1 .*(ic1 <= kmax);
Fhat           = [];

res.IC1            = IC1; 
res.ic1            = ic1; 
res.Fhat           = Fhat0(:, 1 : ic1) * sqrt(T);
res.lambda         = X' * res.Fhat/T;
res.chat           = res.Fhat * res.lambda';

function pos=minindc(x)
% MININDC returns, for each column of x, the row index of its minimum
ncols=size(x,2);
nrows=size(x,1);
pos=zeros(ncols,1);
seq=seqa(1,1,nrows);
for i=1:ncols;
dum=min(x(:,i));
dum1= seq .* ( (x(:,i)-dum) ==0);
pos(i)=sum(dum1);
end;

function x=seqa(a,b,c)
% SEQA returns the column vector a, a+b, ..., a+(c-1)*b
x=linspace(a,(a+(c-1)*b),c)';
