function [fh,xh,gh,H,itct,fcount,retcodeh] = csminwel(fcn,x0,H0,grad,crit,nit,varargin)
%[fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel(fcn,x0,H0,grad,crit,nit,varargin)
% fcn:   string naming the objective function to be minimized
% x0:    initial value of the parameter vector
% H0:    initial value for the inverse Hessian.  Must be positive definite.
% grad:  Either a string naming a function that calculates the gradient, or the null matrix.
%        If it's null, the program calculates a numerical gradient.  In this case fcn must
%        be written so that it can take a matrix argument and produce a row vector of values.
% crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
%        function value by more than crit.
% nit:   Maximum number of iterations.
% varargin: A list of optional length of additional parameters that get handed off to fcn each
%        time it is called.
%        Note that if the program ends abnormally, it is possible to retrieve the current x,
%        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
%        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
%        write g2.mat and g3.mat as well.  If all were written at about the same time, any of them
%        may be a decent starting point.  One can also start from the one with best function value.)
[nx,no]=size(x0);
nx=max(nx,no);
Verbose=1;
NumGrad= isempty(grad);
done=0;
itct=0;
fcount=0;
snit=100;
%tailstr = ')';
%stailstr = [];
% Lines below make the number of Pi's optional.  This is inefficient, though, and precludes
% use of the matlab compiler.  Without them, we use feval and the number of Pi's must be
% changed with the editor for each application.  Places where this is required are marked
% with ARGLIST comments
%for i=nargin-6:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%   stailstr=[' P' num2str(i) stailstr];
%end
f0 = feval(fcn,x0,varargin{:});
%ARGLIST
%f0 = feval(fcn,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
% disp('first fcn in csminwel.m ----------------') % Jinill on 9/5/95
if f0 > 1e50, disp('Bad initial parameter.'), return, end
if NumGrad
   if length(grad)==0
      [g badg] = numgrad(fcn,x0, varargin{:});
      %ARGLIST
      %[g badg] = numgrad(fcn,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
   else
      badg=any(find(grad==0));
      g=grad;
   end
   %numgrad(fcn,x0,P1,P2,P3,P4);
else
   [g badg] = feval(grad,x0,varargin{:});
   %ARGLIST
   %[g badg] = feval(grad,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
end
retcode3=101;
x=x0;
f=f0;
H=H0;
cliff=0;
while ~done
   g1=[]; g2=[]; g3=[];
   %addition fj. 7/6/94 for control
   disp('-----------------')
   disp('-----------------')
   %disp('f and x at the beginning of new iteration')
   disp(sprintf('f at the beginning of new iteration, %20.10f',f))
   %-----------Comment out this line if the x vector is long----------------
      disp([sprintf('x = ') sprintf('%15.8g %15.8g %15.8g %15.8g\n',x)]);
   %-------------------------
   itct=itct+1;
   [f1 x1 fc retcode1] = csminit(fcn,x,f,g,badg,H,varargin{:});
   %ARGLIST
   %[f1 x1 fc retcode1] = csminit(fcn,x,f,g,badg,H,P1,P2,P3,P4,P5,P6,P7,...
   %           P8,P9,P10,P11,P12,P13);
   % itct=itct+1;
   fcount = fcount+fc;
   % erased on 8/4/94
   % if (retcode == 1) | (abs(f1-f) < crit)
   %    done=1;
   % end
   % if itct > nit
   %    done = 1;
   %    retcode = -retcode;
   % end
   if retcode1 ~= 1
      if retcode1==2 | retcode1==4
         wall1=1; badg1=1;
      else
         if NumGrad
            [g1 badg1] = numgrad(fcn, x1,varargin{:});
            %ARGLIST
            %[g1 badg1] = numgrad(fcn, x1,P1,P2,P3,P4,P5,P6,P7,P8,P9,...
            %                P10,P11,P12,P13);
         else
            [g1 badg1] = feval(grad,x1,varargin{:});
            %ARGLIST
            %[g1 badg1] = feval(grad, x1,P1,P2,P3,P4,P5,P6,P7,P8,P9,...
            %                P10,P11,P12,P13);
         end
         wall1=badg1;
         % g1
         save g1 g1 x1 f1 varargin;
         %ARGLIST
         %save g1 g1 x1 f1 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
      end
      if wall1 & (length(H) > 1)% 
         % Bad gradient or back and forth on step length.  Possibly at
         % cliff edge.  Try perturbing search direction if problem not 1D
         %
         %fcliff=fh;xcliff=xh;
         Hcliff=H+diag(diag(H).*rand(nx,1));
         disp('Cliff.  Perturbing search direction.')
         [f2 x2 fc retcode2] = csminit(fcn,x,f,g,badg,Hcliff,varargin{:});
         %ARGLIST
         %[f2 x2 fc retcode2] = csminit(fcn,x,f,g,badg,Hcliff,P1,P2,P3,P4,...
         %     P5,P6,P7,P8,P9,P10,P11,P12,P13);
         fcount = fcount+fc; % put by Jinill
         if  f2 < f
            if retcode2==2 | retcode2==4
                  wall2=1; badg2=1;
            else
               if NumGrad
                  [g2 badg2] = numgrad(fcn, x2,varargin{:});
                  %ARGLIST
                  %[g2 badg2] = numgrad(fcn, x2,P1,P2,P3,P4,P5,P6,P7,P8,...
                  %      P9,P10,P11,P12,P13);
               else
                  [g2 badg2] = feval(grad,x2,varargin{:});
                  %ARGLIST
                  %[g2 badg2] = feval(grad,x2,P1,P2,P3,P4,P5,P6,P7,P8,...
                  %      P9,P10,P11,P12,P13);
               end
               wall2=badg2;
               % g2
               badg2
               save g2 g2 x2 f2 varargin
               %ARGLIST
               %save g2 g2 x2 f2 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
            end
            if wall2
               disp('Cliff again.  Try traversing')
               if norm(x2-x1) < 1e-13
                  f3=f; x3=x; badg3=1;retcode3=101;
               else
                  gcliff=((f2-f1)/((norm(x2-x1))^2))*(x2-x1);
                  if(size(x0,2)>1), gcliff=gcliff', end
                  [f3 x3 fc retcode3] = csminit(fcn,x,f,gcliff,0,eye(nx),varargin{:});
                  %ARGLIST
                  %[f3 x3 fc retcode3] = csminit(fcn,x,f,gcliff,0,eye(nx),P1,P2,P3,...
                  %         P4,P5,P6,P7,P8,...
                  %      P9,P10,P11,P12,P13);
                  fcount = fcount+fc; % put by Jinill
                  if retcode3==2 | retcode3==4
                     wall3=1; badg3=1;
                  else
                     if NumGrad
                        [g3 badg3] = numgrad(fcn, x3,varargin{:});
                        %ARGLIST
                        %[g3 badg3] = numgrad(fcn, x3,P1,P2,P3,P4,P5,P6,P7,P8,...
                        %                        P9,P10,P11,P12,P13);
                     else
                        [g3 badg3] = feval(grad,x3,varargin{:});
                        %ARGLIST
                        %[g3 badg3] = feval(grad,x3,P1,P2,P3,P4,P5,P6,P7,P8,...
                        %                         P9,P10,P11,P12,P13);
                     end
                     wall3=badg3;
                     % g3
                     badg3
                     save g3 g3 x3 f3 varargin;
                     %ARGLIST
                     %save g3 g3 x3 f3 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
                  end
               end
            else
               f3=f; x3=x; badg3=1; retcode3=101;
            end
         else
            f3=f; x3=x; badg3=1;retcode3=101;
         end
      else
         % normal iteration, no walls, or else 1D, or else we're finished here.
         f2=f; f3=f; badg2=1; badg3=1; retcode2=101; retcode3=101;
      end
   else 
      f2=f;f3=f;f1=f;retcode2=retcode1;retcode3=retcode1;
   end
   %how to pick gh and xh
   if f3 < f - crit & badg3==0
      ih=3
      fh=f3;xh=x3;gh=g3;badgh=badg3;retcodeh=retcode3;
   elseif f2 < f - crit & badg2==0
      ih=2
      fh=f2;xh=x2;gh=g2;badgh=badg2;retcodeh=retcode2;
   elseif f1 < f - crit & badg1==0
      ih=1
      fh=f1;xh=x1;gh=g1;badgh=badg1;retcodeh=retcode1;
   else
      [fh,ih] = min([f1,f2,f3]);
      disp(sprintf('ih = %d',ih))
      %eval(['xh=x' num2str(ih) ';'])
      switch ih
         case 1
            xh=x1;
         case 2
            xh=x2;
         case 3
            xh=x3;
      end %case
      %eval(['gh=g' num2str(ih) ';'])
      %eval(['retcodeh=retcode' num2str(ih) ';'])
      retcodei=[retcode1,retcode2,retcode3];
      retcodeh=retcodei(ih);
      if exist('gh')
         nogh=isempty(gh);
      else
         nogh=1;
      end
      if nogh
         if NumGrad
            [gh badgh] = numgrad(fcn,xh,varargin{:});
         else
            [gh badgh] = feval(grad, xh,varargin{:});
         end
      end
      badgh=1;
   end
   %end of picking
   %ih
   %fh
   %xh
   %gh
   %badgh
   stuck = (abs(fh-f) < crit);
   if (~badg)&(~badgh)&(~stuck)
      H = bfgsi(H,gh-g,xh-x);
   end
   if Verbose
      disp('----')
      disp(sprintf('Improvement on iteration %d = %18.9f',itct,f-fh))
   end
   % if Verbose
      if itct > nit
         disp('iteration count termination')
         done = 1;
      elseif stuck
         disp('improvement < crit termination')
         done = 1;
      end
      rc=retcodeh;
      if rc == 1
         disp('zero gradient')
      elseif rc == 6
         disp('smallest step still improving too slow, reversed gradient')
      elseif rc == 5
         disp('largest step still improving too fast')
      elseif (rc == 4) | (rc==2)
         disp('back and forth on step length never finished')
      elseif rc == 3
         disp('smallest step still improving too slow')
      elseif rc == 7
         disp('warning: possible inaccuracy in H matrix')
      end
   % end
   f=fh;
   x=xh;
   g=gh;
   badg=badgh;
end
% what about making an m-file of 10 lines including numgrad.m
% since it appears three times in csminwel.m
