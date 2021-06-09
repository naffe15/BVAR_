function [fX] = bkfilter(X,pl,pu)
%
%	MATLAB COMMAND FOR BAND PASS FILTER: fX = bkfilter(X,pl,pu)
%
%
%	Required Inputs:
%	X     - series of data (T x 1)
%	pl    - minimum period of oscillation of desired component 
%	pu    - maximum period of oscillation of desired component (2<=pl<pu<infinity).
%   
%  Output:
%	fX - matrix (T x 1) containing filtered data 
%   
%   Examples: 
% 		Quarterly data: pl=6, pu=32 returns component with periods between 1.5 and 8 yrs.
%		Monthly data:   pl=2, pu=24 returns component with all periods less than 2 yrs.
%
%	Note:  When feasible, we recommend dropping 2 years of data from the beginning 
%		       and end of the filtered data series.  These data are relatively poorly
%		       estimated.
%
%		===============================================================================
%		This program contains only the default filter recommended in Christiano and 
%		Fitzgerald (1999). This program was written by Eduard Pelz and any errors are 
%		my own and not the authors. For those who wish to use optimal band-pass filters
%		other than the default filter please use the Matlab version of code available 
%		at www.clev.frb.org/Research/workpaper/1999/index.htm next to working paper 9906.
%		===============================================================================
%
%		Version Date: 2/11/00 (Please notify Eduard Pelz at eduard.pelz@clev.frb.org or 
%		(216)579-2063 if you encounter bugs).  
%

if pu <= pl
   error(' (bpass): pu must be larger than pl')
end
if pl < 2
   warning(' (bpass): pl less than 2 , reset to 2')
   pl = 2;
end

[T,nvars] = size(X);

%  This section removes the drift from a time series using the 
%	formula: drift = (X(T) - X(1)) / (T-1).                     
%
undrift = 1;
j = 1:T;
if undrift == 1
   drift = (X(T,1)-X(1,1))/(T-1);
   Xun = X-(j'-1)*drift;
else
   Xun = X;
end

%Create the ideal B's then construct the AA matrix

b=2*pi/pl;
a=2*pi/pu;
bnot = (b-a)/pi;
bhat = bnot/2;

B = ((sin(j*b)-sin(j*a))./(j*pi))';
B(2:T,1) = B(1:T-1,1);
B(1,1) = bnot;

AA = zeros(2*T,2*T);

for i=1:T
   AA(i,i:i+T-1) = B';
   AA(i:i+T-1,i) = B;   
end
AA = AA(1:T,1:T);

AA(1,1) = bhat;
AA(T,T) = bhat;

for i=1:T-1
   AA(i+1,1) = AA(i,1)-B(i,1);
	AA(T-i,T) = AA(i,1)-B(i,1);
end

%Filter data using AA matrix

fX = AA*Xun;