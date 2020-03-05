%% check one: dynamicOutlines


% Goal: Draw tracked ellipses over movie xy31 from 2018-02-01.
%       Label pre-divisional cells in blue and newborn cells in orange.
%       Print particle number, length and width next to tracked cell.



% Strategy:
%
%     0. initialize tracking and image data
%     1. isolate ellipse data from movie (stage xy) of interest
%     2. identify tracks in rejects. all others are tracked.
%     3. for each image, initialize current image
%            4. define lengths, widths, centroids, angles
%            5. draw ellipses from image, colored based on: 
%                 pre-divisional cell = DeepSkyBlue (frame before birth)
%                        newborn cell = Tomato
%           6. display and save
%     7. woohoo!


% last edit: jen, 2019 Mar 5
% commit: edit comments to reflect correct dataset (experiment date)


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


% 1. initialize image data
conversionFactor = 6.5/60;      % scope5 Andor COSMOS = 6.5um pixels / 60x magnification
img_prefix = strcat('lb-fluc-',date,'_xy', num2str(xy), 'T'); 
img_suffix = 'XY1C1.tif';


% 2. load measured data
experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
cd(experimentFolder)
filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
%filename = strcat('lb-fluc-',date,'-width1p7-jiggle-0p5.mat');
load(filename,'D5','T');

    
% 3. compile experiment data matrix
xy_start = xy;
xy_end = xy;
xyData = buildDM(D5, T, xy_start, xy_end, index, expType);
clear D5 T xy_start xy_end e


% 4. trim data to include only full cell cycles
curveFinder = getGrowthParameter(xyData,'curveFinder');   % curve finder (ID of curve in condition)
xyData_fullOnly = xyData(curveFinder > 0,:);
clear curveFinder


% 5. isolate create isDivision vector
isDrop = getGrowthParameter(xyData_fullOnly,'isDrop');   % birth events
div_temp = isDrop(2:end);
isDivision = [div_temp; length(xyData_fullOnly)];
clear isDrop div_temp


% 6. open folder for images of interest (one xy position of experiment)
cd(experimentFolder)
if xy >= 10
    img_folder=strcat('xy', num2str(xy));
else
    img_folder=strcat('xy0', num2str(xy));
end
cd(img_folder);


% 7. create directory of image names in chronological order
imgDirectory = dir(strcat('lb-fluc-',date,'_xy*.tif'));
%imgDirectory = dir(strcat('singleupshift-',date,'_xy*.tif'));
names = {imgDirectory.name};
clear img_folder img_prefix img_suffix experiment newFolder img_folder



% 8. identify tracks and birth events present in each frame
totalFrames = length(imgDirectory); % total frame number
trackIDs = [];
for fr = 1:totalFrames
    
    %tracksInCurrentFrame = xyData_fullOnly(xyData_fullOnly(:,18) == fr,:);    % col 18 = frame #
    tracksInCurrentFrame = xyData_fullOnly(xyData_fullOnly(:,16) == fr,:);                % col 16 = frame #
    trackIDs{fr,1} = tracksInCurrentFrame(:,1);                         % col 1 = TrackID, as given by tracking script ND2Proc_XY
    birthEvents{fr,1} = tracksInCurrentFrame(:,4);                       % col 4 = isDrop
    divEvents{fr,1} = isDivision(xyData_fullOnly(:,16) == fr,:);
end
clear fr tracksInCurrentFrame


% 9. overlay colored cell outlines over each image file
for img = 1:length(names)
    
    % i. initialize current image
    cla
    I=imread(names{img});
    %filename = strcat('dynamicOutlines-widths-fullOnly-xy',num2str(xy),'-frame',num2str(img),'-n',num2str(n),'.tif');
    filename = strcat('dynamicOutlines-birthEvents-xy',num2str(xy),'-frame',num2str(img),'.tif');
    
    figure(1)
    % imtool(I), displays image in grayscale with range
    imshow(I,'DisplayRange',[2200 4800]); %lowering right # increases num sat'd pxls
    
    
    % ii. determine number of cells at birth and pre-division in current frame
    fr_births = birthEvents{img};
    fr_predivs = divEvents{img};
    cells_birth = sum(fr_births);
    cells_prediv = sum(fr_predivs);
    
    
    % iii. skip ellipses if no cells in either category
    if isempty(trackIDs{img}) == 1
        saveas(gcf,filename)
        continue
        
    elseif cells_birth == 0 && cells_prediv == 0  % ii. also, if no birth events 
        saveas(gcf,filename)
        continue
        
    else
        
        % iii. else, isolate data for each image
        dm_currentImage = xyData_fullOnly(xyData_fullOnly(:,16) == img,:);    % col 16 = frame #
        
        majorAxes = getGrowthParameter(dm_currentImage,'length');         %  lengths
        minorAxes = getGrowthParameter(dm_currentImage,'width');          %  widths
        
        centroid_X = getGrowthParameter(dm_currentImage,'x_pos');          % x coordinate of centroid
        centroid_Y = getGrowthParameter(dm_currentImage,'y_pos');          % y coordinate of centroid
        angles = getGrowthParameter(dm_currentImage,'angle');             % angle of rotation of fit ellipses

        IDs = getGrowthParameter(dm_currentImage,'trackID');              % track ID as assigned in ND2Proc_XY
        
        
        % iv. for each cell at birth in current image, draw ellipse in Tomato
        if cells_birth > 0
            
            birth_IDs = IDs(fr_births == 1);
            for pp = 1:length(birth_IDs)
                
                currentP_b = find(IDs == birth_IDs(pp));
                
                [x_rotated, y_rotated] = drawEllipse(currentP_b,majorAxes, minorAxes, centroid_X, centroid_Y, angles, conversionFactor);
                lineVal = 0.5;
                
                color = rgb('Tomato');
                color_text = rgb('DarkRed');
                
                hold on
                plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
                text((centroid_X(currentP_b)+2)/conversionFactor, (centroid_Y(currentP_b)+2)/conversionFactor, strcat('particle:',num2str(IDs(currentP_b))),'Color',color_text,'FontSize',12);
                text((centroid_X(currentP_b)+2)/conversionFactor, (centroid_Y(currentP_b)+6)/conversionFactor, strcat('L:',num2str(majorAxes(currentP_b))),'Color',color_text,'FontSize',12);
                text((centroid_X(currentP_b)+2)/conversionFactor, (centroid_Y(currentP_b)+10)/conversionFactor, strcat('W:',num2str(minorAxes(currentP_b))),'Color',color_text,'FontSize',12);
                xlim([0 2048]);
                ylim([0 2048]);
                
            end
            
        end
        clear pp color birth_IDs currentP_b
   
        
        % v. for each pre-divisional cell at birth in current image, draw ellipse in DeepSkyBlue
        if cells_prediv > 0
            
            prediv_IDs = IDs(fr_predivs == 1);
            for ppp = 1:length(prediv_IDs)
                
                currentP_pd = find(IDs == prediv_IDs(ppp));
                
                [x_rotated, y_rotated] = drawEllipse(currentP_pd,majorAxes, minorAxes, centroid_X, centroid_Y, angles, conversionFactor);
                lineVal = 0.5;
                
                color = rgb('DeepSkyBlue');
                color_text = rgb('Navy');
                
                hold on
                plot(x_rotated,y_rotated,'Color',color,'lineWidth',lineVal)
                text((centroid_X(currentP_pd)+2)/conversionFactor, (centroid_Y(currentP_pd)+2)/conversionFactor, strcat('particle:',num2str(IDs(currentP_pd))),'Color',color_text,'FontSize',12);
                text((centroid_X(currentP_pd)+2)/conversionFactor, (centroid_Y(currentP_pd)+6)/conversionFactor, strcat('L:',num2str(majorAxes(currentP_pd))),'Color',color_text,'FontSize',12);
                text((centroid_X(currentP_pd)+2)/conversionFactor, (centroid_Y(currentP_pd)+10)/conversionFactor, strcat('W:',num2str(minorAxes(currentP_pd))),'Color',color_text,'FontSize',12);
                xlim([0 2048]);
                ylim([0 2048]);
                
            end
            
        end
        
    end
    title(num2str(img))
    
    % 6. save
    saveas(gcf,filename)
    
end



