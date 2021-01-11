function quer(str); 
% function quer(str); 
% Query user 
% if str=='c' Continue or not 
% if str=='g' close Graphs or not
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

