%% check six: average volume accumulation rate vs. birth volume


% Goal: Scatter plot of single-cell average volume accumulation rate
%       vs volume at birth as measured from xy 31 2018-02-01
%
%       ave volume accumulation rate = (Vdiv - Vbirth)/ interdivision time
%       as defined by Dieter. units = cubic um/s


% Strategy:
%
%     0. initialize tracking and image data
%     1. compile experiment data
%     2. trim data to include only full cell cycles
%     3. create division event vector
%     4. isolate volume and interdivision time (tau) data
%     5. isolate data to sizes at birth
%     6. remove data associated with interdivision times of zero
%     7. calculate average volume accumulation rate
%     8. plot scatter and fit line on linear axes (output plot check6-A)
%     9. plot data on log-log scale (output plot check6-B)



% last edit: jen, 2020 April 4
% commit: scatter plot of ave volume accumulation rate vs birth volume


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


% 4. isolate volume and interdivision time (tau) data
volumes = getGrowthParameter(xyData_fullOnly,'volume'); 
cellcycletime = getGrowthParameter(xyData_fullOnly,'curveDurations'); % seconds


% 5. isolate data to sizes at birth and division
birthVol = volumes(birthEvent == 1);
tau = cellcycletime(birthEvent == 1);
divVol = volumes(divEvent ~= 0);


% 6. remove data associated with interdivision times of zero
Vb_true = birthVol(tau > 0);
tau_true = tau(tau > 0);
Vdiv_true = divVol(tau > 0);


% 7. calculate average volume accumulation rate
V_acc = (Vdiv_true-Vb_true)./tau_true;


% 8a. plot scatter
figure(1)
scatter(Vb_true,V_acc)
ylabel('avg vol accumulation rate (um^3/s)')
xlabel('volume at birth (um^3)')
title('single-cell, xy31, 2018-02-01')
axis([2 8 0.001 0.01])


% 8b. fit line to data
fit = polyfit(Vb_true,V_acc,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
m = fit(1);
b = fit(2);
x = linspace(2,8,100);
y = m*x + b;
figure(1)
hold on
plot(x,y)
text(2.5,0.009,strcat('y=',num2str(m),'x+',num2str(b)))


% 9. plot log-log plot of data
figure(2)
scatter(Vb_true,V_acc)
ylabel('avg vol accumulation rate (um^3/s)')
xlabel('volume at birth (um^3)')
title('single-cell, xy31, 2018-02-01')
axis([2 8 0.001 0.01])
set(gca,'xscale','log')
set(gca,'yscale','log')

p = polyfit(log(Vb_true),log(V_acc),1); 
m_loglog = p(1);
b_loglog = exp(p(2));
y_loglog = b_loglog*x.^m_loglog;
figure(2)
hold on
plot(x,y_loglog)
text(2.2,0.009,'y=b*x^m')
text(2.2,0.0082,strcat('b=',num2str(b_loglog)))
text(2.2,0.0076,strcat('m=',num2str(m_loglog)))
