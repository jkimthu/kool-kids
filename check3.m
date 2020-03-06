%% check three: cell size distributions


% Goal: Plot scatter of single-cell width vs volume at birth
%       as tracked from xy31 of 2018-02-01.


% Strategy:
%
%     0. initialize tracking and image data
%     1. compile experiment data
%     2. trim data to include only full cell cycles
%     3. isolate width, volume and birth event data
%     4. isolate data to sizes at birth
%     5. plot scatter!



% last edit: jen, 2020 Mar 5
% commit: scatter plot of single-cell width vs volume


% OK LET'S GO!


%% Initialize experiment data

clc
clear

% 0. initialize complete meta data
cd('/Users/jen/Documents/StockerLab/Data_analysis/')
load('storedMetaData.mat')

% 0. initialize experiment and xy movie to analyze
index = 15; % 2018-02-01
xy = 31;

% 0. initialize experiment meta data
date = storedMetaData{index}.date;
expType = storedMetaData{index}.experimentType;


% 0. load measured data
experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
cd(experimentFolder)
filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
load(filename,'D5','T');

    
% 1. compile experiment data
xy_start = xy;
xy_end = xy;
xyData = buildDM(D5, T, xy_start, xy_end, index, expType);
clear D5 T xy_start xy_end e


% 2. trim data to include only full cell cycles
curveFinder = getGrowthParameter(xyData,'curveFinder');   % curve finder (ID of curve in condition)
xyData_fullOnly = xyData(curveFinder > 0,:);
clear curveFinder


% 3. isolate width, volume and birth event data  
widths = getGrowthParameter(xyData_fullOnly,'width'); 
volumes = getGrowthParameter(xyData_fullOnly,'volume'); 
birthEvent = getGrowthParameter(xyData_fullOnly,'isDrop'); 


% 4. isolate data to sizes at birth
birthWid = widths(birthEvent == 1);
birthVol = volumes(birthEvent == 1);


% 5. plot scatter
figure(1)
scatter(birthWid,birthVol)
ylabel('Width at birth (um)')
xlabel('Volume at birth (cubic um)')
title('Single-cell width vs vol, xy31, 2018-02-01')
xlim([1 1.8])





