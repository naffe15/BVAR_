function [ir,Omeg] = iresponse_zeros_signs(Phi,Sigma,hor,lag,var_pos,f,sr,draws,toler)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'iresponse_zeros_sign' computes the impulse response functions using zero
% and sign restrictions on the endogenous variables
% References: 
% Arias, J. E., Rubio-Ram?rez, J. F. and Waggoner, D. F.: 2018, Inference
% Based on SVARs Identified with Sign and Zero Restrictions: Theory and
% Applications, Econometrica 86, 685–720. 
% Binning, A.: 2013, Underidentified SVAR models: A framework for combining
% short and long-run restrictions with sign-restrictions, Working Paper
% 2013/14, Norges Bank.  

% Inputs:
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - hor, horizon of the IRF
% - unit, 1 shock STD or 1 percent increase
% - (var_pos,f,sr) inputs for the zero and sign restrictions see bvar.m and
% tutorial_.m 

% Output:
% - ir contains the IRF 
% 1st dimension:   variable 
% 2st dimension:   horizon 
% 3st dimension:   shock

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin< 8
    draws = 1; % Number of draws
end
if nargin< 9
    toler = 10000; % Number of rotation attempt
end

p  = lag;
k  = size(Sigma,1);
C1 = chol(Sigma,'lower');

ir      = nan(k,hor,k);
Omeg    = nan(k);

T = hor; % length of impulse response function. 

shocks = eye(k);

[Q,index,flag] = findQs(k,f);

if flag == 1
    error('Rank condition not satisfied, the model is overidentified');
end

shock_pos = logical(shocks); % position of shock

% var_pos = [1,1,4,1,1];
% var_pos = [1,2,2,3,3];       % position of corresponding variable to shock 
% % eg monetary policy shock should result in a positive increase in interest 
% % rates, aggregate demand shock should result in a positive increase in gdp
% % etc.

R = zeros(k,T,length(shocks),draws); % Contains impulse resonse functions
% 1st dimension = variable
% 2nd dimension = time
% 3th dimension = shock
% 4rd dimension = draw

counter = 1;

B      = [Phi(end,:); Phi(1:end-1,:)];
Btilde = B(2:end,:)';
% Btilde = Phi(1:end-1,:)';
alpha = [Btilde;eye(k*(p-1)),zeros(k*(p-1),k)]; % Build companion form matrix
if draws > 10
    wb = waitbar(0, 'Generating Rotations');
end

tj = 0;
while counter < draws+1
    
    tj = tj +1;
    
    C = generateDraw(C1,k);
    
    P = findP(C,B,Q,p,k,index);
    
    W = C*P;

    for jj = 1:length(shocks)
        
        if W(var_pos(jj),jj) < 0
            shock = -shocks(:,jj);
        else
            shock = shocks(:,jj);
        end
        
        V = zeros(k*p,T);
        
        V(1:k,1) = W*shock;
        
        chk = W*shock;
        sr_index = ~isnan(sr(:,jj));
        tmp = sign(chk(sr_index)) - sr(sr_index,jj);

        if any(tmp~=0)
            jj = 0;
            break
        end
        
        for ii = 2:T
            V(:,ii) = alpha*V(:,ii-1);
        end
        
        R(:,:,jj,counter) = V(1:k,:);
        
    end
    
    if jj == length(shocks)
        counter = counter + 1;
    end
    
    if tj > toler
        warning('I could not find a rotation')
        return;
    end
    if draws > 10, waitbar(counter/draws, wb); end
end


if draws > 10, close(wb); end

Omeg = W; 
ir    = R;