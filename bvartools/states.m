function [st] = states(bcp5,bct5,nbp,nbt,nd)
%STATES Expansion and contraction indicator from peaks and troughs.
%
%   st = STATES(bcp5,bct5,nbp,nbt,nd)
%
%   Inputs
%     bcp5, bct5   positions of the peaks and of the troughs (from RAWALL)
%     nbp, nbt     number of peaks and of troughs
%     nd           number of observations
%
%   Output
%     st           nd x 1   1 expansion, 0 contraction
%
%   See also DATE_, RAWALL.

% 1=expansion;0=contraction  

st=zeros(nd,1);

if bcp5(1)<bct5(1)

	i=2;
	% need to determine if pattern ends with [trough;peak]  assume ends with peak 
	
	nrp=min([nbt;nbp])';
	nrpp=nrp+1;
	
	if nbt==nbp
		%  ends with trough
		nrpp=nrp;
    end
	
	
	st(1:bcp5(1))=ones(bcp5(1),1);
	while i<=nrpp
     
		st(bct5(i-1)+1:bcp5(i))=ones(bcp5(i)-bct5(i-1),1);
		i=i+1;
    end
	
	
	% test if no peaks equal no troughs
	% pattern then is pt,pt,pt etc & fill in from trough to end with ones 
	
	
	if nbt==nrpp
		st(bct5(nrp)+1:nd)=ones(nd-bct5(nrp),1);
    end

else

	i=1;
	nrp=nbp;
	while i<=nrp
	
	%  bct5[i];
	%  bcp5[i];
	
	st(bct5(i)+1:bcp5(i))=ones(bcp5(i)-bct5(i),1);
	i=i+1;
    end
	
	
	% test if one more troughs than peaks if so then pattern is tptpt & need to fill in with ones at end
	
	
	if nbt>nbp
      	st(bct5(nbp+1)+1:nd)=ones(nd-bct5(nbp+1),1);
    end


end



% End of function
