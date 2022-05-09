
%% Fetch Data
%Fetch data sets up a function so that Haver Data can be imported directly
%to Matlab

%Inputs:
% Database - 
% seriesName - used for the acutal names coming from Haver
% startDate and endDate are previously chosen for the range of data you'd like to pull

function [data,timev]=fetchData(database,seriesName,startDate,endDate,frequency)
        D=fetch(database,seriesName,startDate,endDate,frequency);
        data=D(:,2);
        timev =D(:,1);
end