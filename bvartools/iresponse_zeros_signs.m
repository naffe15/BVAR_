function [ir,Omeg] = iresponse_zeros_signs(Phi,Sigma,hor,lag,var_pos,f,sr,draws,toler,normalization_)
%IRESPONSE_ZEROS_SIGNS Responses with zero (short and long run) and impact sign restrictions.
%
%   [ir,Omeg] = IRESPONSE_ZEROS_SIGNS(Phi,Sigma,hor,lag,var_pos,f,sr)
%   [ir,Omeg] = IRESPONSE_ZEROS_SIGNS(...,draws,toler,normalization_)
%
%   Uses the algorithm of Binning (2013): each candidate rotation satisfies
%   the zero restrictions by construction and is kept if the impact signs
%   are right. Used by BVAR_ with options.zeros_signs.
%
%   Inputs
%     Phi             (ny*lag+nx) x ny  VAR coefficients, lag-1 block first
%     Sigma           ny x ny   covariance of the reduced-form residuals
%     hor             scalar    number of horizons
%     lag             scalar    number of lags
%     var_pos         1 x ny    for each shock, the variable whose impact
%                               response is normalised to be positive
%                               (used only if normalization_ = 1)
%     f, sr           restriction matrices from SIGN2MATRIX
%     draws           scalar    accepted rotations to return (default 1)
%     toler           scalar    max number of candidates (default 10000)
%     normalization_  0 or 1    sign normalisation by var_pos (default 0)
%
%   Outputs
%     ir     ny x hor x ny x draws   responses (variable x horizon x shock
%                                    x draw); NaN if no rotation is found
%     Omeg   ny x ny                 last accepted rotation of
%                                    chol(Sigma,'lower')
%
%   The importance weight of Arias, Rubio-Ramirez and Waggoner (2018) is
%   not applied. Errors if the zero restrictions over-identify the model.
%
%   See also SIGN2MATRIX, FINDQS, FINDP, GENERATEDRAW, IRESPONSE_SET.
%
%   References: Arias, J. E., J. F. Rubio-Ramirez and D. F. Waggoner (2018),
%   Econometrica 86(2), 685-720. Binning, A. (2013), "Underidentified SVAR
%   models: A framework for combining short and long-run restrictions with
%   sign-restrictions", Norges Bank Working Paper 2013/14.

%   Filippo Ferroni, 6/1/2015.

if nargin < 8
    draws = 1; % Number of draws
end
if nargin < 9
    toler = 10000; % Number of rotation attempt
end
if nargin < 10
    normalization_ = 0; % sign normalization on rotations
end
p         = lag;
k         = size(Sigma,1);
endo_numb = p*k;
C1        = chol(Sigma,'lower');

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

B      = [Phi(endo_numb+1:end,:); Phi(1:endo_numb,:)];
Btilde = Phi(1:endo_numb,:)';
% B      = [Phi(end,:); Phi(1:end-1,:)];
% Btilde = B(2:end,:)';
% Btilde = Phi(1:end-1,:)';
alpha = [Btilde;eye(k*(p-1)),zeros(k*(p-1),k)]; % Build companion form matrix
WW    = nan(k);
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
                
        shock = shocks(:,jj);
        if normalization_ == 1
            if W(var_pos(jj),jj) < 0
                shock = -shocks(:,jj);
            end
        end
        
        WW(:,jj) = W*shock;

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

%Omeg = W; 
Omeg = C1\WW; 
ir    = R;