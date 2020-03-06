%% check four: single cell volume at division vs birth


% Goal: Scatter plot of single-cell volume at division vs volume at birth.


% Strategy:
%
%     0. initialize tracking and image data
%     1. compile experiment data
%     2. trim data to include only full cell cycles
%     3. create division event vector
%     4. isolate length, width, volume and birth event data
%     5. isolate data to sizes at birth
%     6. plot scatter!
%     7. fit line to data



% last edit: jen, 2020 Mar 5
% commit: scatter plot of single-cell div vs birth volume


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


% 3. create division event vector
birthEvent = getGrowthParameter(xyData_fullOnly,'isDrop');   % birth events
div_temp = birthEvent(2:end);
divEvent = [div_temp; length(xyData_fullOnly)];
clear div_temp


% 4. isolate length, width, volume and birth event data
lengths = getGrowthParameter(xyData_fullOnly,'length');   
widths = getGrowthParameter(xyData_fullOnly,'width'); 
volumes = getGrowthParameter(xyData_fullOnly,'volume'); 


% 5. isolate data to sizes at birth and division
birthLen = lengths(birthEvent == 1);
birthWid = widths(birthEvent == 1);
birthVol = volumes(birthEvent == 1);

divLen = lengths(divEvent ~= 0);
divWid = widths(divEvent ~= 0);
divVol = volumes(divEvent ~= 0);


% 5. plot scatter
figure(1)
scatter(birthVol,divVol)
ylabel('volume at division (cubic um)')
xlabel('volume at birth (cubic um)')
title('single-cell, xy31, 2018-02-01')


% 6. fit line to data
fit = polyfit(birthVol,divVol,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
m = fit(1);
b = fit(2);
x = linspace(0,14,100);
y = m*x + b;
figure(1)
hold on
plot(x,y)
text(10,20,strcat('y=',num2str(m),'x+',num2str(b)))

