
%% Fetch Data
%Fetch data sets up a function so that Haver Data can be imported directly
%to Matlab

%Inputs:
% Database - refers to Haver database such as USECON, USNA, or CAPSTOCK; USNA
% most commonly used to keep consistency with Jeff's previous chain aggregation code
% and other BEA numbers.
% seriesName - used for the acutal names coming from Haver
% startDate and endDate are previously chosen for the range of data you'd like to pull

function [data,timev]=fetchData(database,seriesName,startDate,endDate,frequency)
        D=fetch(database,seriesName,startDate,endDate,frequency);
        data=D(:,2);
        timev =D(:,1);
end