function [x,rc] = csolve(FUN,x,gradfun,crit,itmax,varargin)
%function [x,rc] = csolve(FUN,x,gradfun,crit,itmax,varargin)
% FUN should be written so that any parametric arguments are packed in to x,
% and so that if presented with a matrix x, it produces a return value of
% same dimension of x.  The number of rows in x and FUN(x) are always the
% same.  The number of columns is the number of different input arguments
% at which FUN is to be evaluated.
%
% gradfun:  string naming the function called to evaluate the gradient matrix.  If this
%           is null (i.e. just "[]"), a numerical gradient is used instead.
% crit:     if the sum of absolute values that FUN returns is less than this,
%           the equation is solved.
% itmax:    the solver stops when this number of iterations is reached, with rc=4
% varargin: in this position the user can place any number of additional arguments, all
%           of which are passed on to FUN and gradfun (when it is non-empty) as a list of 
%           arguments following x.
% rc:       0 means normal solution, 1 and 3 mean no solution despite extremely fine adjustments
%           in step length (very likely a numerical problem, or a discontinuity). 4 means itmax
%           termination.
%---------- delta --------------------
% differencing interval for numerical gradient
delta = 1e-6;
%-------------------------------------
%------------ alpha ------------------
% tolerance on rate of descent
alpha=1e-3;
%---------------------------------------
%------------ verbose ----------------
verbose=1;% if this is set to zero, all screen output is suppressed
%-------------------------------------
%------------ analyticg --------------
analyticg=1-isempty(gradfun); %if the grad argument is [], numerical derivatives are used.
%-------------------------------------
nv=length(x);
tvec=delta*eye(nv);
done=0;
if isempty(varargin)
   f0=feval(FUN,x);
else
   f0=feval(FUN,x,varargin{:});
end   
af0=sum(abs(f0));
af00=af0;
itct=0;
while ~done
   if itct>3 & af00-af0<crit*max(1,af0) & rem(itct,2)==1
      randomize=1;
   else
      if ~analyticg
         if isempty(varargin)
            grad = (feval(FUN,x*ones(1,nv)+tvec)-f0*ones(1,nv))/delta;
         else
            grad = (feval(FUN,x*ones(1,nv)+tvec,varargin{:})-f0*ones(1,nv))/delta;
         end
      else % use analytic gradient
         grad=feval(gradfun,x,varargin{:});
      end
      if isreal(grad)
         if rcond(grad)<1e-12
            grad=grad+tvec;
         end
         dx0=-grad\f0;
         randomize=0;
      else
         if(verbose),disp('gradient imaginary'),end
         randomize=1;
      end
   end
   if randomize
      if(verbose),fprintf(1,'\n Random Search'),end
      dx0=norm(x)./randn(size(x));
   end
   lambda=1;
   lambdamin=1;
   fmin=f0;
   xmin=x;
   afmin=af0;
   dxSize=norm(dx0);
   factor=.6;
   shrink=1;
   subDone=0;
   while ~subDone
      dx=lambda*dx0;
      f=feval(FUN,x+dx,varargin{:});
      af=sum(abs(f));
      if af<afmin
         afmin=af;
         fmin=f;
         lambdamin=lambda;
         xmin=x+dx;
      end
      if ((lambda >0) & (af0-af < alpha*lambda*af0)) | ((lambda<0) & (af0-af < 0) )
         if ~shrink
            factor=factor^.6;
            shrink=1;
         end
         if abs(lambda*(1-factor))*dxSize > .1*delta;
            lambda = factor*lambda;
         elseif (lambda > 0) & (factor==.6) %i.e., we've only been shrinking
            lambda=-.3;
         else %
            subDone=1;
            if lambda > 0
               if factor==.6
                  rc = 2;
               else
                  rc = 1;
               end
            else
               rc=3;
            end
         end
      elseif (lambda >0) & (af-af0 > (1-alpha)*lambda*af0)
         if shrink
            factor=factor^.6;
            shrink=0;
         end
         lambda=lambda/factor;
      else % good value found
         subDone=1;
         rc=0;
      end
   end % while ~subDone
   itct=itct+1;
   if(verbose)
      fprintf(1,'\nitct %d, af %g, lambda %g, rc %g',itct,afmin,lambdamin,rc)
      fprintf(1,'\n   x  %10g %10g %10g %10g',xmin);
      fprintf(1,'\n   f  %10g %10g %10g %10g',fmin);
   end
   x=xmin;
   f0=fmin;
   af00=af0;
   af0=afmin;
   if itct >= itmax
      done=1;
      rc=4;
   elseif af0<crit;
      done=1;
      rc=0;
   end
end
