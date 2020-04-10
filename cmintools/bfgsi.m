function H = bfgsi(H0,dg,dx)
% H = bfgsi(H0,dg,dx)
% dg is previous change in gradient; dx is previous change in x;
% 6/8/93 version that updates inverse hessian instead of hessian
% itself.
% Copyright by Christopher Sims 1996.  This material may be freely
% reproduced and modified.
if size(dg,2)>1
   dg=dg';
end
if size(dx,2)>1
   dx=dx';
end
Hdg = H0*dg;
dgdx = dg'*dx;
if (abs(dgdx) >1e-12)
   H = H0 + (1+(dg'*Hdg)/dgdx)*(dx*dx')/dgdx - (dx*Hdg'+Hdg*dx')/dgdx;
else
   disp('bfgs update failed.')
   disp(['|dg| = ' num2str(sqrt(dg'*dg)) '|dx| = ' num2str(sqrt(dx'*dx))]);
   disp(['dg''*dx = ' num2str(dgdx)])
   disp(['|H*dg| = ' num2str(Hdg'*Hdg)])
   H=H0;
end
save H.dat H
