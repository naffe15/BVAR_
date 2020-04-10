struct csminwel_functions_parameters {
  gsl_matrix * H;
  parameters * p;
}
  function [fh,xh,gh,H,itct,fcount,retcodeh] = csminwel(gsl_multimin_function_fdf * fcn,x0,H0,grad,,varargin)
