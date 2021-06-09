function [sindex] = OrderingIndexes(varordering,varnames,newstrng)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'OrderingIndexes' find the specified order of variables

% input: 
% - varnames    = the list (names) of all the variables in the database
% - varordering = the names and the order of the variables in the VAR

% output: 
% - sindex = the index in the varnames that correspond to the var ordering
% - index_XX = the index that the variable XXX has in the VAR

% Filippo Ferroni, 
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    newstrng='';
end
    
[~,sindex] = ismember(varordering,varnames);
for ii = 1 : size(varordering,2) 
    check=0;
    for jjj = 1: size(varnames,2)        
         if strcmp(deblank(varnames{jjj}),deblank(varordering{ii})) == 1,
             check = 1;
         end             
    end
    if check == 1
        assignin('base', [newstrng 'index_' varordering{ii} ], ii);
    else
        warning(['I did not find ' varordering{ii} ' in varnames']) 
    end
end
