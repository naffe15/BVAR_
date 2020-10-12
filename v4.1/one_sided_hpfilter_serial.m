function [ytrend,ycycle]=one_sided_hpfilter_serial(y,lambda,discard)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% I am grateful to Martin Kliem for pointing out this alternative to the
% Kalman variant to me and for sharing his code with me.
%
% Copyright: Alexander Meyer-Gohde
%
% You are free to use/modify/redistribute this program so long as original
% authorship credit is given and you in no way impinge on its free
% distribution
%
%
% This program executes a one-sided HP filter by running the
% standard two-sided HP filter successively with each new observation
%
%  Input: y - a Txn data matrix, where T is the number of observations
%            on n variables (i.e.,data is assumed to be in column format). 
%         lambda - a smoothing scalar. If not entered, default is 1600.
%         discard  -  a scalar. The first discard periods will be
%            discarded resulting in output matrices of size
%            (T-discard)xn. Optional: if not entered, default value is 0.
%
% Output: ytrend - a (T-discard)xn matrix of extracted trends for
%                  each of the n variables.
%         ycycle a  (T-discard)xn matrix of deviations from the trends 
%                     for  each of the n variables. Optional.
%
%       Usage example:
%       [ytrend]=one_sided_hp_filter_kalman(y)
%                   will yield the Txn matrix of trends using the data in
%                   y with lambda set to 1600
%
%
% This one-side HP filter finds a series {ytrend_t}_{t=1}^T for
% each n by calculating for each t the standard HP filtered trend using
% all data upto that t and setting ytrend_t equal to the resulting trend 
% value for period t.
%
% See: Mehra, Y.P. (2004). "The Output Gap, Expected Future Inflation and 
%                        Inflation Dynamics: Another Look," 
%           The B.E. Journal of Macroeconomics, Berkeley Electronic Press.
%
% The standard HP filter finds a series {ytrend_t}_{t=1}^T that solves
% the following minimization probelm
% 
% min sum_{t=1}^T(y_t-ytrend_t)
% +lambda*sum_{t=2}^{T-1}[(ytrend_{t+1}-ytrend_{t})-(ytrend_{t}-ytrend_{t-1})]
%
% Rearanging the first order conditions yields:  A*ytrend=y   where
% A=[   1+lambda    -2*lambda       lambda        0               0           0 ...
%     [    -2*lambda   1+5*lambda   -4*lambda       lambda       0           0 ...
%     [    lambda       -4*lambda        1+6*lambda   -4*lambda  lambda  0 ...
%     [                                         ...
%     [    0      ...   0   lambda       -4*lambda        1+6*lambda   -4*lambda  lambda 0  ... 0]
%     [                                         ...      ]
%     [  0   ...    lambda       -4*lambda        1+6*lambda  -4*lambda  lambda]
%     [  0   ...           0     lambda  -4*lambda       1+5*lambda      -2*lambda ]
%     [  0   ...           0       0        lambda              -2*lambda    1+lambda]
%
%
%   Hodrick, R.J. and E.C.Prescott (1997), "Postwar U.S. Business Cycles:
%   An Empirical Investigation." Jounal of Money, Credit and Banking.
%   29(1), Feb. pp. 1--16.
%
%  Note that filter will be calculated over and over again with an
%  expanding matrix A. This program uses sparse matrices and exploits the
%  patter in A as t progresses.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,  lambda = 1600; end%If the user does not provide lambda, set it to 1600
[T,n] = size (y); % Calculate the number of periods and the number of variables in the series

%%Preliminary calculations
x1=[1+lambda, -2*lambda, lambda]; %The non-zero elements of the first row of A
x2=[ -2*lambda, 1+5*lambda, -4*lambda, lambda];%The non-zero elements of the second row of A
x3=[lambda, -4*lambda, 1+6*lambda, -4*lambda, lambda];%The non-zero elements of thej'th row of A, 2<j<T-2
x2rev=x2(end:-1:1);%The non-zero elements of the second-to-last row of A (just x2 in reverse)
x1rev=x1(end:-1:1); %The non-zero elements of the last row of A (just x1 in reverse)

%The minimization problem for the HP filter is either trivially given by
%the following or is not defined when there are less than three
%observations
for t=1:2
    ytrend(t,:)=y(t,:); %Thus, the trend is set to the observation
end
t=t+1;
%For three observations, the minimization problem yields the following
junk=[x1;-2*lambda, 1+4*lambda, -2*lambda;x1rev]\y(1:t,:);%Build the matrix A and solve the system A*y_standard_hp_trend=y
ytrend(t,:)=junk(t,:);%select the t'th element of the resulting trend
%Starting with four observations, the pattern using the preliminary
%calculations begins
t=t+1;
Ibegin=[1; 1; 1; 2; 2; 2; 2];%make a list (a column vector) containing at position i the row of the i'th non-zero element in the first two rows of A 
Jbegin=[1;2;3;1;2;3;4];%make a list (a column vector) containing at position i the column of the i'th non-zero element in the first two rows of A 
Xbegin=[x1';x2'];%make a list (a column vector) containing at position i the i'th non-zero element in the first two rows of A 
Iend=[3; 3; 3; 3; 4; 4; 4];%make a list (a column vector) containing at position i the row of the i'th non-zero element in the last two rows of A 
Jend=[1;2;3;4;2;3;4];%make a list (a column vector) containing at position i the column of the i'th non-zero element in the last two rows of A 
Xend=[x2rev';x1rev'];%make a list (a column vector) containing at position i the i'th non-zero element in the last two rows of A 
junk=sparse([Ibegin; Iend],[Jbegin; Jend],[Xbegin; Xend])\y(1:t,:);
ytrend(t,:)=junk(t,:);
%For the fifth observation on up, the simple pattern can be exploited by
%adding the new row and column indicies of the new non-zero elements to the Ibegin and Jbegin as well as
%the values of the elements themselves. Note that this requires the
%row and column indicies of the non-zeor elements in the last two rows to
%be advanced
while t<T
    t=t+1;
    Ibegin=[Ibegin;t-2;t-2;t-2;t-2;t-2];%Add to the list the rows of the non-zero elements of the t-2'th row
    Jbegin=[Jbegin;t-4;t-3;t-2;t-1;t];%Add to the list the col.umns of the non-zero elements of the t-2'th row
    Xbegin=[Xbegin;x3'];%Add to the list the non-zero elements of the t-2'th row
    Iend=Iend+1;%Advance the rows of the non-zero elements in the last two rows
    Jend=Jend+1;%Advance the columns of the non-zero elements in the last two rows
    junk=sparse([Ibegin; Iend],[Jbegin; Jend],[Xbegin; Xend],t,t,(t-4)*5+14)\y(1:t,:);%Build the matrix A and solve the system A*y_standard_hp_trend=y
    ytrend(t,:)=junk(t,:);%select the t'th element of the resulting trend
end

if nargout==2%Should the user have requested a second output
    ycycle=y-ytrend;%The second output will be the deviations from the HP trend
end

if nargin==3%If the user provided a discard parameter
    ytrend=ytrend(discard+1:end,:);%Remove the first "discard" periods from the trend series
    if nargout==2%Should the user have requested a second output
        ycycle=ycycle(discard+1:end,:);%The second output will be the deviations from the HP trend, likewise with  first "discard" periods removed
    end
end
end