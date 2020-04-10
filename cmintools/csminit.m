function [fhat,xhat,fcount,retcode] = csminit(fcn,x0,f0,g0,badg,H0,varargin)
% [fhat,xhat,fcount,retcode] = csminit(fcn,x0,f0,g0,badg,H0,...
%                                       P1,P2,P3,P4,P5,P6,P7,P8)
% retcodes: 0, normal step.  5, largest step still improves too fast.
% 4,2 back and forth adjustment of stepsize didn't finish.  3, smallest
% stepsize still improves too slow.  6, no improvement found.  1, zero
% gradient.
%---------------------
% Modified 7/22/96 to omit variable-length P list, for efficiency and compilation.
% Places where the number of P's need to be altered or the code could be returned to
% its old form are marked with ARGLIST comments.
%
% Fixed 7/17/93 to use inverse-hessian instead of hessian itself in bfgs
% update.
%
% Fixed 7/19/93 to flip eigenvalues of H to get better performance when
% it's not psd.
%
%tailstr = ')';
%for i=nargin-6:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%end
%ANGLE = .03;
ANGLE = .005;
%THETA = .03;
THETA = .3; %(0<THETA<.5) THETA near .5 makes long line searches, possibly fewer iterations.
FCHANGE = 1000;
MINLAMB = 1e-9;
% fixed 7/15/94
% MINDX = .0001;
% MINDX = 1e-6;
MINDFAC = .01;
fcount=0;
lambda=1;
xhat=x0;
f=f0;
fhat=f0;
g = g0;
gnorm = norm(g);
%
if (gnorm < 1.e-12) & ~badg % put ~badg 8/4/94
   retcode =1;
   dxnorm=0;
   % gradient convergence
else
   % with badg true, we don't try to match rate of improvement to directional
   % derivative.  We're satisfied just to get some improvement in f.
   %
   %if(badg)
   %   dx = -g*FCHANGE/(gnorm*gnorm);
   %  dxnorm = norm(dx);
   %  if dxnorm > 1e12
   %     disp('Bad, small gradient problem.')
   %     dx = dx*FCHANGE/dxnorm;
   %   end
   %else
   % Gauss-Newton step;
   %---------- Start of 7/19/93 mod ---------------
   %[v d] = eig(H0);
   %toc
   %d=max(1e-10,abs(diag(d)));
   %d=abs(diag(d));
   %dx = -(v.*(ones(size(v,1),1)*d'))*(v'*g);
%      toc
   dx = -H0*g;
%      toc
   dxnorm = norm(dx);
   if dxnorm > 1e12
      disp('Near-singular H problem.')
      dx = dx*FCHANGE/dxnorm;
   end
   dfhat = dx'*g0;
   %end
   %
   %
   if ~badg
      % test for alignment of dx with gradient and fix if necessary
      a = -dfhat/(gnorm*dxnorm);
      if a<ANGLE
         dx = dx - (ANGLE*dxnorm/gnorm+dfhat/(gnorm*gnorm))*g;
         % suggested alternate code:  ---------------------
         dx = dx*dxnorm/norm(dx)    % This keeps scale invariant to the angle correction
         % ------------------------------------------------
         dfhat = dx'*g;
         % dxnorm = norm(dx);  % this line unnecessary with modification that keeps scale invariant
         disp(sprintf('Correct for low angle: %g',a))
      end
   end
   disp(sprintf('Predicted improvement: %18.9f',-dfhat/2))
   %
   % Have OK dx, now adjust length of step (lambda) until min and
   % max improvement rate criteria are met.
   done=0;
   factor=3;
   shrink=1;
   lambdaMin=0;
   lambdaMax=inf;
   lambdaPeak=0;
   fPeak=f0;
   lambdahat=0;
   while ~done
      if size(x0,2)>1
         dxtest=x0+dx'*lambda;
      else
         dxtest=x0+dx*lambda;
      end
      % home
      f = feval(fcn,dxtest,varargin{:});
      %ARGLIST
      %f = feval(fcn,dxtest,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
      % f = feval(fcn,x0+dx*lambda,P1,P2,P3,P4,P5,P6,P7,P8);
      disp(sprintf('lambda = %10.5g; f = %20.7f',lambda,f ))
      %debug
      %disp(sprintf('Improvement too great? f0-f: %g, criterion: %g',f0-f,-(1-THETA)*dfhat*lambda))
      if f<fhat
         fhat=f;
         xhat=dxtest;
         lambdahat = lambda;
      end
      fcount=fcount+1;
      shrinkSignal = (~badg & (f0-f < max([-THETA*dfhat*lambda 0]))) | (badg & (f0-f) < 0) ;
      growSignal = ~badg & ( (lambda > 0)  &  (f0-f > -(1-THETA)*dfhat*lambda) );
      if  shrinkSignal  &   ( (lambda>lambdaPeak) | (lambda<0) )
         if (lambda>0) & ((~shrink) | (lambda/factor <= lambdaPeak))
            shrink=1;
            factor=factor^.6;
            while lambda/factor <= lambdaPeak
               factor=factor^.6;
            end
            %if (abs(lambda)*(factor-1)*dxnorm < MINDX) | (abs(lambda)*(factor-1) < MINLAMB)
            if abs(factor-1)<MINDFAC
               if abs(lambda)<4
                  retcode=2;
               else
                  retcode=7;
               end
               done=1;
            end
         end
         if (lambda<lambdaMax) & (lambda>lambdaPeak)
            lambdaMax=lambda;
         end
         lambda=lambda/factor;
         if abs(lambda) < MINLAMB
            if (lambda > 0) & (f0 <= fhat)
               % try going against gradient, which may be inaccurate
               lambda = -lambda*factor^6
            else
               if lambda < 0
                  retcode = 6;
               else
                  retcode = 3;
               end
               done = 1;
            end
         end
      elseif  (growSignal & lambda>0) |  (shrinkSignal & ((lambda <= lambdaPeak) & (lambda>0)))
         if shrink
            shrink=0;
            factor = factor^.6;
            %if ( abs(lambda)*(factor-1)*dxnorm< MINDX ) | ( abs(lambda)*(factor-1)< MINLAMB)
            if abs(factor-1)<MINDFAC
               if abs(lambda)<4
                  retcode=4;
               else
                  retcode=7;
               end
               done=1;
            end
         end
         if ( f<fPeak ) & (lambda>0)
            fPeak=f;
            lambdaPeak=lambda;
            if lambdaMax<=lambdaPeak
               lambdaMax=lambdaPeak*factor*factor;
            end
         end
         lambda=lambda*factor;
         if abs(lambda) > 1e20;
            retcode = 5;
            done =1;
         end
      else
         done=1;
         if factor < 1.2
            retcode=7;
         else
            retcode=0;
         end
      end
   end
end
disp(sprintf('Norm of dx %10.5g', dxnorm))
