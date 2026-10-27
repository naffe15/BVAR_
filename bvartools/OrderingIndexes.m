function [sindex] = OrderingIndexes(varordering,varnames,newstrng)
%ORDERINGINDEXES Positions of selected variables and workspace index variables.
%
%   sindex = ORDERINGINDEXES(varordering,varnames)
%   sindex = ORDERINGINDEXES(varordering,varnames,newstrng)
%
%   Inputs
%     varordering   cell   names of the VAR variables, in the VAR order
%     varnames      cell   names of all variables in the data set
%     newstrng      char   prefix of the index variables (default '')
%
%   Output
%     sindex        position of each element of varordering in varnames
%
%   Also creates in the base workspace, for each variable XXX, a variable
%   [newstrng 'iXXX'] equal to its position in the VAR (spaces and
%   hyphens removed from the name).

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
        %assignin('base', [newstrng 'index_' varordering{ii} ], ii);
        String = varordering{ii};
        String = replace(String,'-','');
        iString = String(find(~isspace(String)));
        assignin('base', [newstrng 'i' iString ], ii);
    else
        warning(['I did not find ' varordering{ii} ' in varnames']) 
    end
end
