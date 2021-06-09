%function [varargout] = date_(varargin)
function     [dt_]  = ...
    date_(x, time, strfreq, t0, t1, options) 
%    date_(x, time, freq, ndfy, ndfm, ndly, ndlm, turnphase, phase, cycle, thresh, nrep, complete)
% st,trinary,notentp,dura_,ampl_,cumm_,excc_,durcv_,amplcv_,exccv_
%***********************************************************************************************/
%*                                                                                             */
%*  Modified BBQ program                                                                       */
%*                                                                                             */
%*  Computes turning points and imposes restrictions in one step                               */
%*                                                                                             */
%*    NB: rule is peak greater than turnphase on either side,                                  */
%*             trough less than turnphase on either side.           (can be modified)          */
%*    Restriction: Phase quarters/months, Cycle  quarters/months Alternating Troughs and Peaks,*/
%*                  if two peaks(troughs) in a row chooses highest(lowest), Trough must be     */
%*                  lower then preceeding peak, and No turning point within phase length       */
%*                  of end points                                                              */
%*                                                                                             */
%*                                                                                              */
%*                                                                                             */
%*    Date: 15th November 2005                                                                  */
%*     Author: James Engel - Code modified from Adrian Pagan and Don Harding's BBQ code     */
%*
%*     further adjustments  made  by  F. Canova
%***********************************************************************************************/
% 

ndfy = t0(1); % year start
ndfm = t0(2); % month/quarter start, 
ndly = t1(1); % year start, 
ndlm = t1(2); % month/quarter start, 

if length(x) ~= length(time)
    error('data (x) and time must have the same length');
end

if strcmp(strfreq,'q') == 1 
    freq =1;
    % default settings cycle parameters q data
    turnphase   = 2;
    phase       = 2;        % censoring rules
    cycle       = 5;        % lenght  of cycle
    
    if ndfm>4 || ndlm>4
        error('Something went wrong. Frequency is quarterly; therefore, starting and ending quarters cannot be larger than 4!');
    end
    
elseif strcmp(strfreq,'m') == 1 
    freq =2;
    % default settings cycle parameters m data
    turnphase   = 6;
    phase       = 6;        % censoring rules
    cycle       = 15;        % lenght  of cycle
    
else
    error('strfreq must be ''q'' (quarterly) or ''m'' (monthly)');
end

notentp=0;

if freq==1 % q
    nd= 4*(ndly-1-ndfy) + (5-ndfm) + (ndlm);     % Number of data points %
    stepsize = 1/4;
elseif freq==2 % m
    nd= 12*(ndly-1-ndfy) + (13-ndfm) + (ndlm);
    stepsize = 1/12;
end

thresh      = 10.4;     % bypasses phase and cycle restriction if peak to trough is > than thresh
nrep        = 1;        % 1 if analyze real data
complete    = 1;        % if= 1- use complete cycles,if =0 -use incomplete cycles (excess still on complete cycle)

if nargin > 5
    if isfield(options,'turnphase')==1
        turnphase = options.turnphase;
    end
    % censoring rules 
    if isfield(options,'phase')==1
        phase = options.phase;
    end
    % lenght of cycle
    if isfield(options,'cycle')==1
        cycle = options.cycle;
    end
    % bypasses phase and cycle restriction if peak to trough is > than thresh
    if isfield(options,'thresh')==1
        thresh = options.thresh;
    end
    % 1 if analyze real data
    if isfield(options,'nrep')==1
        nrep = options.nrep;
    end
    % if= 1- use complete cycles,if =0 -use incomplete cycles (excess still on complete cycle)
    if isfield(options,'complete')==1
        complete = options.complete;
    end
    
end

pdcv=0;tdcv=0;pacv=0;tacv=0;pdm=0;pdma=0;tdm=0;tdma=0;pdcm=0;tdcm=0;
pdem=0;tdem=0;tdema=0;pdema=0;epcv=0;etcv=0;

[bcp5, bct5, nbp, nbt] = rawall(x(1:nd),turnphase,nd,phase,cycle,thresh);   % calculates turning points with restrictions %

if nbp+nbt<=2
    notentp=notentp+1;
else
    
    %ntr=nbt;npk=nbp;
    
    nr  = [nbt;nbp];
    nv  = max(nr);
    
    %               pdc=zeros(nv,1);tdc=zeros(nv,1);
    %               pda=zeros(nv,1);tda=zeros(nv,1);td=zeros(nv,1);pd=zeros(nv,1);
    pdc = zeros(nv,1);
    tdc = zeros(nv,1); % pde=zeros(nv,1);tde=zeros(nv,1);
    %               tdea=zeros(nv,1);pdea=zeros(nv,1);
    
    % calculate peak to trough durations & amps
    % p st&s for peaks,t for troughs,p gives cntractions,t expansions
    % code is pd,td is durations; pda,tda is amps; pdc,tdc is cum move pde, tde exces
    % excess is measured differently to avoid case that amps is close to zero in part cycle so that denom can become neg.so use cum movements as denom
    % now vs triangle in early paper
    
    if bcp5(1,1) < bct5(1,1)      % Peaks are first
        
        nr  = [nbt;nbp];
        r   = nbt;
        pd  = bct5(1:r,1)-bcp5(1:r,1);
        pda = x(bct5(1:r,1))-x(bcp5(1:r,1));
        k   = 1;
        
        while k<=r
            %pdc(k)=sumc(x(bcp5(k,1):bct5(k,1),1)-x(bcp5(k,1),1));
            pdc(k)=sum(x(bcp5(k,1):bct5(k,1),1)-x(bcp5(k,1),1));
            k=k+1;
        end
    else                      % troughs are First
        
        r   = nbt-1;
        pd  = bct5(2:r+1,1)-bcp5(1:r,1);
        pda = x(bct5(2:r+1,1))-x(bcp5(1:r,1));
        
        k=1;
        while k<=r
            %pdc(k)=sumc(x(bcp5(k,1):bct5(k+1,1),1)-x(bcp5(k,1),1));
            pdc(k)=sum(x(bcp5(k,1):bct5(k+1,1),1)-x(bcp5(k,1),1));
            k=k+1;
        end
        
        %r1 = r;
        
    end
    
    % calculate trough to peak durations & amplitudes
    
    if bct5(1,1) < bcp5(1,1)        %  Troughs are first
        r   = nbp;
        td  = bcp5(1:r,1)-bct5(1:r,1);
        tda = x(bcp5(1:r,1))-x(bct5(1:r,1));
        k   = 1;
        while k<=r
            %tdc(k)=sumc(x(bct5(k,1):bcp5(k,1),1)-x(bct5(k,1),1));
            tdc(k)=sum(x(bct5(k,1):bcp5(k,1),1)-x(bct5(k,1),1));
            k=k+1;
        end
        
    else                      % peaks are First
        r   = nbp-1;
        td  = bcp5(2:r+1,1)-bct5(1:r,1);
        tda = x(bcp5(2:r+1,1))-x(bct5(1:r,1));
        
        
        k=1;
        while k<=r
            
            %tdc(k)=sumc(x(bct5(k,1):bcp5(k+1,1),1)-x(bct5(k,1),1));
            tdc(k)=sum(x(bct5(k,1):bcp5(k+1,1),1)-x(bct5(k,1),1));
            k=k+1;
        end
        
    end
    pdc=pdc(1:size(pd,1));
    tdc=tdc(1:size(td,1));
    
    % compute excesses
    %excess is percentage of triangle area
    za      = (pd.*pda)/2;
    pde     = 100*(pdc-za-.5*pda)./za;
    
    pdea    = 100*(pdc-((pd.*pda)/2))./za;    
    za      = (td.*tda)/2;
    tde     = 100*(tdc-za-.5*tda)./za;
    
    tdea    = 100*(tdc-((td.*tda)/2))./za;
    
    %***********************************************************************************************/
    
    
    
    if complete == 0    % switch: 1-only use complete cycles,0-use incomplete cycles (excess still on complete cycle)%
        
        bct5u   = bct5;
        bcp5u   = bcp5;
        
        if bcp5(1,1) < bct5(1,1)     % modifies code to include incomplete cycles %
            
            bct5u  = [1;bct5];
            
        elseif 	bct5(1,1) < bcp5(1,1)
            
            bcp5u = [1;bcp5];
            
        end
        
        
        nbtu = size(bct5u,1);
        nbpu = size(bcp5u,1);
        
        if bcp5u(nbpu,1) < bct5u(nbtu,1)    % modifies code to include incomplete cycles %
            
            bcp5u  = [bcp5u;nd];
            
        elseif 	bct5u(nbt,1) < bcp5u(nbp,1)
            
            bct5u = [bct5u;nd];
            
        end
        
        
        ntr     = size(bct5u,1);
        npk     = size(bcp5u,1);
        % nr      = [ntr;npk];
        % nv=max(nr);
        
        
        %        pdc=zeros(nv,1);tdc=zeros(nv,1);
        %        pda=zeros(nv,1);tda=zeros(nv,1);td=zeros(nv,1);pd=zeros(nv,1);
        %        pdc=zeros(nv,1);tdc=zeros(nv,1);
        
        
        if bcp5u(1,1) < bct5u(1,1)       % Peaks are first
            
            %nr  = [ntr;npk];
            r   = ntr;
            pd  = bct5u(1:r,1)-bcp5u(1:r,1);
            pda = x(bct5u(1:r,1))-x(bcp5u(1:r,1));
            k   = 1;
            while k<=r
                %pdc(k)=sumc(x(bcp5u(k,1):bct5u(k,1),1)-x(bcp5u(k,1),1));
                pdc(k)=sum(x(bcp5u(k,1):bct5u(k,1),1)-x(bcp5u(k,1),1));
                k=k+1;
            end
        else                      % troughs are First
            r   = ntr-1;
            pd  = bct5u(2:r+1,1)-bcp5u(1:r,1);
            pda = x(bct5u(2:r+1,1))-x(bcp5u(1:r,1));
            k   = 1;
            while k<=r
                pdc(k)=sum(x(bcp5u(k,1):bct5u(k+1,1),1)-x(bcp5u(k,1),1));
                %pdc(k)=sumc(x(bcp5u(k,1):bct5u(k+1,1),1)-x(bcp5u(k,1),1));
                k=k+1;
            end
            
%             r1=r;
            
        end
        
        % calculate trough to peak durations & amplitudes
        
        if bct5u(1,1) < bcp5u(1,1)       %  Troughs are first
            
            r   = npk;
            td  = bcp5u(1:r,1)-bct5u(1:r,1);
            tda = x(bcp5u(1:r,1))-x(bct5u(1:r,1));
            k   = 1;
            while k<=r
                %tdc(k)=sumc(x(bct5u(k,1):bcp5u(k,1),1)-x(bct5u(k,1),1));
                tdc(k)=sum(x(bct5u(k,1):bcp5u(k,1),1)-x(bct5u(k,1),1));
                k=k+1;
            end
            
            
        else                     % peaks are First
            r   = npk-1;
            td  = bcp5u(2:r+1,1)-bct5u(1:r,1);
            tda = x(bcp5u(2:r+1,1))-x(bct5u(1:r,1));
            
            
            k=1;
            while k<=r                
                %tdc(k)=sumc(x(bct5u(k,1):bcp5u(k+1,1),1)-x(bct5u(k,1),1));
                tdc(k)=sum(x(bct5u(k,1):bcp5u(k+1,1),1)-x(bct5u(k,1),1));                
                k=k+1;
            end
            
        end
        pdc=pdc(1:size(pd,1));
        tdc=tdc(1:size(td,1));
        
    end
    
    
    
    %***********************************************************************************************/
    
    
    % cumulate.when nrep=1 then gives raw data,otherwsise sums over monte carlo amounts
    
    pdcv    = pdcv + stdc(pd)/meanc(pd);       % compute cv's
    tdcv    = tdcv + stdc(td)/meanc(td);
    pacv    = pacv + stdc(pda)/meanc(pda);
    tacv    = tacv + stdc(tda)/meanc(tda);
    epcv    = epcv + stdc(pde)/meanc(pde);  %  cv of xss
    etcv    = etcv + stdc(tde)/meanc(tde);
    
    pdm     = pdm+meanc(pd);          % durations
    pdma    = pdma+meanc(pda);
    
    tdm     = tdm+meanc(td);          % durations
    tdma    = tdma+meanc(tda);
    
    pdcm    = pdcm+meanc(pdc);       % cumulative
    tdcm    = tdcm+meanc(tdc);
    
    % pdem=pdem+meanc(pde);
    tdem    = tdem+meanc(tde);       % xss
    tdema   = tdema+meanc(tdea);
    pdem    = pdem+meanc(pde);
    pdema   = pdema+meanc(pdea);
   
end

if nrep==1&&(nbp+nbt>2)
    
    %nbp and nbt=number of peaks and number of troughs %
    
    disp('peaks at')
    fprintf('%.2f\n', time(bcp5))
    %disp(bcp5(1:nbp))
    
    disp('troughs at')
    fprintf('%.2f\n', time(bct5))
    %disp(bct5(1:nbt))
    
    % total  number of  turning  points
    ttp=length(bcp5)+length(bct5);
    
    % collecting  turning  points in  trinary  indicator
    trinary=zeros(length(x),1);
    trinary(bcp5)=1;
    trinary(bct5)=-1;
    
    
    
end


if nrep==1&&(nbp+nbt>2)
    %remove below if want states ...
    %st=zeros(nd,1);
    
    [st] = states(bcp5,bct5,nbp,nbt,nd);
    dt_.st(:,nrep) = st;
    
    %determine obs in which states have been completed%...
    na  = min([bct5(1);bcp5(1)])';
    nb  = max([bct5(nbt);bcp5(nbp)])';
    %nb-na+1;
    z   = [x(na:nb) st(na:nb)];
    
    %disp('dated series')
    %z1= [time(time_start+1:time_end)' x st];
    %format  short  g
    %round(z1, 2)
    %%fprintf('%.2f\n', '%.2f\n', '%.2f\n', z1)
    
end
nrep1   = nrep-notentp;

% Compute final statistics
% duration contractions/duration expansions
dt_.dura = [pdm/nrep1 tdm/nrep1];
% amplitudes contractions/amplitude expansions
dt_.ampl  = 100*[pdma/nrep1 tdma/nrep1];
% cumulative contractions/cumulative expansions
dt_.cumm  = 100*[pdcm/nrep1 tdcm/nrep1];
% excess movements percent of triangle area
% contractions/expansions
dt_.excc  = [pdem/nrep1 tdem/nrep1];
% cv of durations contractions/expansions
dt_.durcv = [pdcv/nrep1 tdcv/nrep1];
% cv of amplitude contractions/expansions
dt_.amplcv = [pacv/nrep1 tacv/nrep1];
% cv of excess movements contractions/expansions
dt_.exccv = [epcv/nrep1 etcv/nrep1];
% 
dt_.trinary = trinary;
% 
dt_.notentp = notentp;


function m=meanc(x)
% get mean of each column
m=mean(x);
if size(m,2)>1
   m=m';
end

function m=stdc(x)
% get standard deviation of each column
m=std(x);
if size(m,2)>1
   m=m';
end
