function [f_hat_save,f_hat_save_separate, grid, step] = pw_cohesion(dY,q_low, q_high,factor,idx)
% Purpose: Computes Power Cohesion
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% SOURCE:
% Financial cycles: Characterisation and Real-Time Measurement by 
% By Yves Schuler, Paul Hiebert, and Tuomas Peltonen (2020)
%--------------------------------------------------------------------------
% Code Originally Written by Yves Schuler.
% Code Edited by Rob Wlodarski
% Supervised by Fabio Canova
% Code Last Updated: 14 September 2023
%--------------------------------------------------------------------------
% INPUT:
% dY - time series to be estimated (N x Y), 
% N= # of observations, Y = # of variables (max=8)
% qlow - lower bound for spectral density (e.g., 5)
% qhigh - upper bound for spectral density (e.g., 200)
% factor - precision of spectral density estimates (e.g., 8)
% idx - indices of the variables of interest (e.g., [1 1 1 1])
%--------------------------------------------------------------------------
% OUTPUT:
% f_hat_save - power cohesion measure.
% f_hat_save_separate - cross-spectra for the chosen variables.
% grid - designed grid (to be useful for graphs & tables)
% step - step (for frequency resolution in graphs).

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
if nargin < 5
    idx = ones(size(dY,2),1);
end

%--------------------------------------------------------------------------
% Part 1: Defining grid & step. 
omega_step = 2048;
resolution = 16;
% Defining grid in different cases. 
if q_low==0 && q_high ~=0
    grid = (2*pi/(q_high)):2*pi/omega_step:pi-pi/(omega_step);
elseif q_low ~=0 && q_high ==0
    grid = 0:2*pi/omega_step:(2*pi)/q_low;
elseif q_low==0 && q_high ==0
    grid = 0:2*pi/omega_step:pi-pi/(omega_step);
else    
    grid = (2*pi/(q_high)):2*pi/omega_step:(2*pi/q_low);
end
% Define steps for frequency display in graphs
step = (max(grid) - min(grid))/resolution;                    

%--------------------------------------------------------------------------
% Part 2: Choosing variables of interest & designing placeholders. 
dy = dY(:,logical(idx));
M = size(dy,2);
f_hat_save = zeros(length(grid),1);
f_hat_save_separate = [];

%--------------------------------------------------------------------------
% Part 3: Computing power cohesion. 
kk = 0;
for ii=1:M-1
    for jj=ii+1:M
    kk = kk + 1;
    y=dy(:,ii);
    x=dy(:,jj);
    inv=(x~=x)|(y~=y); % clean missing or invalid data points
    x(inv)=[];
    y(inv)=[];        
    T = length(y);
    gammahat = xcov(y,x)/T;  
    gammahat = gammahat./(std(y)*std(x));
    t_cov = length(gammahat);
    Mwind = round(sqrt(length(y))*factor);                     
    w = parzenwin(2*Mwind+1);
    ti=floor(t_cov/2);
    tt=(ti+Mwind):-1:(ti-Mwind);
    cov_trunc = gammahat(tt,:).*(w);
    f_hat=(2/(2*pi))*abs((cov_trunc'*exp(-1i*(-Mwind:Mwind)'*grid))'); 
    f_hat_save = f_hat_save + f_hat;
    f_hat_save_separate = [f_hat_save_separate f_hat.*(std(y)*std(x))];
    end
end
f_hat_save = f_hat_save/kk;
end

