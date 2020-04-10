% initialize the MH by minimizing the - logposterior kernel
T = length(y);
vSwitch = ones(2, 1);
% % The parameters might be bounded, but the minimization routine will not
% % care, thus we need to make sure that the routine can go anywhere and
% % that we are still able to stay within the bounds.
% x0 = boundsINV(guess);
x0 = guess;
% posterior maximization
[fh, xh, gh, H, itct, fcount, retcodeh] = csminwel('logpostkernel',x0,.01*eye(length(x0)),[],10e-5,1000,centered,y, vSwitch, iPriorScale);
% processing the output of the maximization
postmode = xh;
JJ       = jacob(xh);
HH       = JJ * H * JJ';