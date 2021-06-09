function [fX,AA] = cffilter(X,pl,pu,root,drift,ifilt,nfix,thet)
%
%	MATLAB COMMAND FOR DEFAULT FILTER:  fX = cffilter(X,pl,pu)
%
%	Required Inputs:
%  	X		- matrix of data (T x # variables) with variables in columns
%		pl		- minimum period of oscillation of desired component 
%		pu		- maximum period of oscillation of desired component (2<=pl<pu<infinity).
%
%  Output:
%		fX 	- matrix (T x nvars) containing filtered data 
%
%	Examples:  
%   Quarterly data: pl=6, pu=32 returns component with periods between 1.5 and 8 yrs.
%	 Monthly data  : pl=2, pu=24 returns component with all periods less than 2 yrs.
%
%  Note:  When feasible, we recommend dropping 2 years of data from the beginning 
%		and end of the filtered data series.  These data are relatively poorly estimated.
%
%	Optional Inputs:  root,drift,ifilt,nfix,thet.  THESE INPUTS ARE NOT REQUIRED.
%
%	Additional Output:  AA - T x T matrix which premultiplies a T x Nvars data matrix 
%				to compute the band-pass filtered data.  (Note: data should have any drift
%				term or time trend removed before applying AA.)
%
%=======================================================================================
%		The next discussion is only for users who wish to use optimal band pass filters 
%		other than the default filter recommended in Christiano and Fitzgerald (1999). 
%=======================================================================================
%
%  Version Date: 6/10/99 (Please notify Terry Fitzgerald at
%			tj.fitzgerald@clev.frb.org or (216)579-2293 if you encounter bugs).  
%
%  This is a matlab function that filters time series data using  
%  approximations to the band pass filter as discussed in the paper
%  "The Band Pass Filter" by Lawrence J. Christiano and Terry J. Fitzgerald (1999).

%	Several band-pass approximation strategies can be selected in this 
%  program.  To explain how to select these alternatives, we rely upon the notation 
%  developed in the paper.  The default setting used above returns the filtered data fX 
%  associated with the unrestricted optimal filter assuming a Random Walk with drift
%  time series process.  
%
%	Alternative filters require selecting the type of filter (asymetric, symmetric, or 
%  fixed length symmetric) and the underlying time series assumption (random walk, IID, 
%  or a user provided MA process with or without a unit root).  The fixed length 
%  symmetric filter advocated by Baxter and King(1999), and the trigonometric regression
%	approach, can also be chosen.  
%
%  Optional Inputs:
%		root		= 0	no unit root in time series						(default: root = 1);
%					= 1	unit root in time series
%		drift		= 0	no drift or time trend in time series			(default: drift = 1);
%					= 1	time series is Random Walk with possibly a nonzero drift
%					 		or time series is stationary about a linear time trend.
%		ifilt 	= 0	Asymmetric Filter 									(default: ifilt = 0);
%					= 1	Symmetric Filter
%					= 2	Fixed Length Symmetric Filter
%					= 3 	Baxter-King Fixed Length Symmetric Filter 
%					= 4	Trigonometric Regression Filter
%		nfix	= sets fixed length (used w/ ifilt=2 or 3)				(default: nfix=-1).
%		thet	= MA coefficients for time series model 					(default: thet=1)
%				if root=1 :	x(t) = mu + x(t-1) + thet(1)e(t) + thet(2)e(t-1) + ...
%				if root=0 :	x(t) = thet(1)e(t) + thet(2)e(t-1) + ...
%						Examples:
%									root=1, thet=1 	Random Walk filter
%									root=0, thet=1		IID filter
%
%  The following function calls are permitted:
%		fX = cffilter(X,pl,pu)										: optimal Random Walk filter
%		fX = cffilter(X,pl,pu,root)								: choose unit root option
%		fX = cffilter(X,pl,pu,root,drift)						: choose drift/trend option
%		fX = cffilter(X,pl,pu,root,drift,ifilt)				: choose type of filter
%		fX = cffilter(X,pl,pu,root,drift,ifilt,nfix)			: choose fixed filter length
%		fX = cffilter(X,pl,pu,root,drift,ifilt,nfix,thet)	: choose MA time series model
%
%  Notes:  
%	1) Users wishing to use a fixed length symmetric filter (ifilt = 2 or 3) must
%     also supply the lead/lag length of the filter (nfix).
%	2) The filtered data matrix (fX)	associated with some filters include zeros at the 
%		beginning and end of the data set.  The filter is not available for these points.
%		You should truncate these points from your analysis.  Examples include
%		fixed length filters and filters for user supplied time series with length(thet)>1.  
%		The default filter returns all nonzero data.
%
if nargin < 4
   root = 1;
   drift = 1;
   ifilt = 0;
   nfix = -1;
   thet = 1;
elseif nargin < 5
   drift = 1;
   ifilt = 0;
   nfix = -1;
   thet = 1;
elseif nargin < 6
   ifilt = 0;
   nfix = -1;
   thet = 1;
elseif nargin < 7
   nfix = -1;
   thet = 1;
elseif nargin < 8
   thet = 1;
end

nq=length(thet)-1;
b=2*pi/pl;
a=2*pi/pu;

[T,nvars] = size(X);

if T < 5
   warning(' (bpass): T < 5')
end
if T < 2*nq+1
   error(' (bpass): T must be at least 2*q+1');
end
if pu <= pl
   error(' (bpass): pu must be larger than pl')
end
if pl < 2
   warning(' (bpass): pl less than 2 , reset to 2')
   pl = 2;
end
if root ~= 0 && root ~= 1
   error(' (bpass): root must be 0 or 1')
end
if drift<0 || drift > 1
   error(' (bpass): drift must be 0 or 1')
end
if ifilt<0 || ifilt > 4
   error(' (bpass): ifilt must be 0, 1, 2, 3, or 3')
end
if (ifilt == 2 || ifilt == 3) && nfix < 1
   error(' (bpass): fixed lag length must be >= 1')
end
if ifilt == 2 && nfix < nq
   error(' (bpass): fixed lag length must be >= q')
end
if (ifilt == 2 || ifilt ==3) && nfix >= T/2
      error(' (bpass): fixed lag length must be < T/2')
end
if (ifilt == 4 && T-2*floor(T/2) ~= 0)
   error(' (bpass): trigonometric regressions only available for even T')
end

[m1,m2]=size(thet);
if m1 > m2 
   th=thet;
else
   th=thet';
end
%   compute g(thet)
%   [g(1) g(2) .... g(2*nq+1)] correspond to [c(q),c(q-1),...,c(1),
%                                        c(0),c(1),...,c(q-1),c(q)]
%   cc = [c(0),c(1),...,c(q)]
thp=flipud(th);
g=conv(th,thp);
cc = g(nq+1:2*nq+1);
%   compute "ideal" Bs
j = 1:2*T;
B = [(b-a)/pi (sin(j*b)-sin(j*a))./(j*pi)]';
%    compute R using closed form integral solution
R = zeros(T,1);
if nq > 0
   R0 = B(1)*cc(1) + 2*B(2:nq+1)'*cc(2:nq+1);
   R(1) = pi*R0;
   for i = 2:T
      dj = Bge(i-2,nq,B,cc);
      R(i) = R(i-1) - dj;
   end
else
   R0 = B(1)*cc(1);
   R(1) = pi*R0;
   for j = 2:T
      dj = 2*pi*B(j-1)*cc(1);
      R(j) = R(j-1) - dj;
   end
end
%
%======================================================================
%
AA = zeros(T,T);
%
%----------------------------------
%  asymmetric filter
%----------------------------------
%
if ifilt == 0
   if nq==0
      for i = 1:T
         AA(i,i:T) = B(1:T-i+1)';
         if root == 1
            AA(i,T) = R(T+1-i)/(2*pi);
         end
      end
      AA(1,1) = AA(T,T);
      %  Use symmetry to construct bottom 'half' of AA
      AA = AA + flipud(fliplr(triu(AA,1)));
   else
      % CONSTRUCT THE A MATRIX size T x T
      A = Abuild(T,nq,g,root);
      Ainv = inv(A);
      % CONSTRUCT THE d MATRIX size T x 1 
      for np = 0:ceil(T/2-1)
         d = zeros(T,1);
         ii = 0;
        	for jj = np-root:-1:np+1+root-T
           	ii = ii+1;
	        d(ii) = Bge(jj,nq,B,cc);
   	        end
         if root == 1
            d(T-1) = R(T-np);
         end
         %  COMPUTE Bhat = inv(A)*d
         Bhat = Ainv*d;
         AA(np+1,:) = Bhat';
      end
      %  Use symmetry to construct bottom 'half' of AA
      AA(ceil(T/2)+1:T,:) = flipud(fliplr(AA(1:floor(T/2),:)));
   end
end   
%
%-------------------------------------
%  symmetric filter
%-------------------------------------
%
if ifilt == 1
   if nq==0
      for i = 2:ceil(T/2)
         np = i-1;
         AA(i,i:i+np) = B(1:1+np)';
         if root == 1
            AA(i,i+np) = R(np+1)/(2*pi);
         end
         AA(i,i-1:-1:i-np) = AA(i,i+1:i+np);
      end
      %  Use symmetry to construct bottom 'half' of AA
      AA(ceil(T/2)+1:T,:) = flipud(fliplr(AA(1:floor(T/2),:)));
   else
      for np = nq:ceil(T/2-1)
         nf = np;
         nn = 2*np+1;
      % CONSTRUCT THE A MATRIX size nn x nn
       	A = Abuild(nn,nq,g,root);
         Ainv = inv(A);
      % CONSTRUCT THE d MATRIX size nn x 1
         d = zeros(nn,1);
         ii = 0;
         for jj = np-root:-1:-nf+root
  	         ii = ii+1;
     	      d(ii) = Bge(jj,nq,B,cc);
         end
         if root == 1
            d(nn-1) = R(nf+1);
	     end            
      %  COMPUTE Bhat = inv(A)*d
         Bhat = Ainv*d;
         AA(np+1,1:2*np+1) = Bhat';
      end
      %  Use symmetry to construct bottom 'half' of AA
      AA(ceil(T/2)+1:T,:) = flipud(fliplr(AA(1:floor(T/2),:)));
   end
end
%
%----------------------------------
%  fixed length symmetric filter
%----------------------------------
%
if ifilt == 2
   if nq==0
      bb = zeros(2*nfix+1,1);
      bb(nfix+1:2*nfix+1) = B(1:nfix+1);
      bb(nfix:-1:1) = B(2:nfix+1);
      if root == 1
         bb(2*nfix+1) = R(nfix+1)/(2*pi);
         bb(1) = R(nfix+1)/(2*pi);
      end
      for i = nfix+1:T-nfix
         AA(i,i-nfix:i+nfix) = bb';
      end
   else
      nn = 2*nfix+1;
      % CONSTRUCT THE A MATRIX size nn x nn
     	A = Abuild(nn,nq,g,root);
      Ainv = inv(A);
      % CONSTRUCT THE d MATRIX size nn x 1
      d = zeros(nn,1);
      ii = 0;
      for jj = nfix-root:-1:-nfix+root
 	      ii = ii+1;
	  	   d(ii) = Bge(jj,nq,B,cc);
      end
      if root == 1
     	   d(nn-1) = R(nn-nfix);
      end
      %  COMPUTE Bhat = inv(A)*d
      Bhat = Ainv*d;
      for ii = nfix+1:T-nfix
         AA(ii,ii-nfix:ii+nfix) = Bhat';
      end
   end
end
%
%----------------------------------
%  Baxter-King filter
%----------------------------------
%
if ifilt == 3
   bb = zeros(2*nfix+1,1);
   bb(nfix+1:2*nfix+1) = B(1:nfix+1);
   bb(nfix:-1:1) = B(2:nfix+1);
   bb = bb - sum(bb)./(2*nfix+1);
   for i = nfix+1:T-nfix
      AA(i,i-nfix:i+nfix) = bb';
   end
end
%
%----------------------------------
%  Trigonometric Regression filter
%----------------------------------
%
if ifilt == 4
	jj = 1:T./2;
	% find frequencies in desired band omitting T/2;
	jj = find(T./pu<=jj & jj<=T./pl & jj<T./2);
	if isempty(jj)
	   error(' (bpass): frequency band is empty in trigonometric regression');
	end
	om = 2*pi*jj./T;
	if pl > 2
  		for t = 1:T
   		for k = T:-1:1
      		l = t-k;
      		tmp = sum(cos(om*l));
      		AA(t,k) = tmp;
   		end
  		end
	else
  		for t = 1:T
   		for k = T:-1:1
      		l = t-k;
      		tmp = sum(cos(om*l));
      		tmp2 = ( cos(pi*(t-l)).*cos(pi*t) )./ 2;
      		AA(t,k) = tmp + tmp2;
   		end
  		end
	end
	AA = AA*2/T;
end
%
%======================================================================
%
%  check that sum of all filters equal 0 if assuming unit root
%
if root == 1
   tst = max(abs(sum(AA')));
   if tst > 1.e-09 && root ~= 0
      warning(' (bpass): Bhat does not sum to 0 ')
      tst   
   end
end
%
%======================================================================
%
%	compute filtered time series using selected filter matrix AA
%
if drift > 0
   X = undrift(X);
end
fX = AA*X;
%
%======================================================================
%				Functions
%======================================================================
%
function y = Bge(jj,nq,B,cc)
%
%  closed form solution for integral of B(z)g(z)(1/z)^j  (eqn 16)
%     nq > 0, jj >= 0
%     if nq = 0, y = 2*pi*B(absj+1)*cc(1);
%
absj =abs(jj);
if absj >= nq
   dj = B(absj+1)*cc(1) + B(absj+2:absj+nq+1)'*cc(2:nq+1);
   dj = dj + flipud(B(absj-nq+1:absj))'*cc(2:nq+1);
elseif absj >= 1
   dj = B(absj+1)*cc(1) + B(absj+2:absj+nq+1)'*cc(2:nq+1);
   dj = dj + flipud(B(1:absj))'*cc(2:absj+1);
   dj = dj + B(2:nq-absj+1)'*cc(absj+2:nq+1);
else
   dj = B(absj+1)*cc(1) + 2*B(2:nq+1)'*cc(2:nq+1);
end
y = 2*pi*dj;
%
%-----------------------------------------------------------------------
%
function A = Abuild(nn,nq,g,root)
%
%  builds the nn x nn A matrix in A.12
%   if root == 1 (unit root)
%      Abig is used to construct all but the last 2 rows of the A matrix
%   elseif root == 0 (no unit root)
%      Abig is used to construct the entire A matrix
%
if root == 1
	Abig=zeros(nn-2,nn+2*(nq-1));
	for j = 1:nn-2      
	   Abig(j,j:j+2*nq) = g'; 
	end
	A = Abig(:,nq:nn+nq-1);
	%   construct A(-f)
	Q = -ones(nn-1,nn);
	Q = tril(Q);
	F = zeros(1,nn-1);
	F(nn-1-nq:nn-1) = g(1:nq+1);
	A(nn-1,:) = F*Q;
	%    construct last row of A
   A(nn,:) = ones(1,nn);
else
	Abig=zeros(nn,nn+2*(nq-1));
	for j = 1:nn      
   	Abig(j,j:j+2*nq) = g'; 
	end
	A = Abig(:,nq+1:nn+nq);
end   
%    multiply A by 2*pi
A = 2*pi*A;
%
%-----------------------------------------------------------------------
%
function xun = undrift(x)
%
%  This function removes the drift or a linear time trend from a time series using the formula
%			drift = (x(T) - x(1)) / (T-1).
%
%  Input:  x - data matrix x where columns represent different variables, x is (T x # variables). 
%  Output: xun - data matrix same size as x with a different drift/trend removed from each variable.
%
[T,nvars] = size(x);
xun = zeros(T,nvars);
dd = (0:T-1)';
for ivar = 1:nvars
   drift = (x(T,ivar)-x(1,ivar)) / (T-1);
   xun(:,ivar) = x(:,ivar) - dd*drift;
end

