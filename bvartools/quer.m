function quer(str); 
%QUER Asks the user, in a dialog box, whether to continue or close figures.
%
%   QUER('c')   asks whether to continue; 'No' stops with an error
%   QUER('g')   asks whether to close all figures
%
%   See also KFILTERNAN.

str=upper(str); 
switch str 
case 'C' 
    
    button=questdlg('Continue? Y/N',...
        'Continue Program',...
        'Yes','No','Yes'); 
    switch button 
    case 'No' 
        error('Program execution stopped'); 
    case 'Yes' 
        fprintf('%-30%s \n','  '); 
    end
case 'G'
    button=questdlg('Close All Graphs Y/N?','Graphs','No','Yes','Yes'); 
    switch button 
    case 'No' 
        fprintf('%-30%s \n','  '); 
    case 'Yes' 
        close all; 
    end
end

