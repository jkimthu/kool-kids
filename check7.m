%% check seven: cropped single cells for Dieter (predivisional)


% Goal: Modify check5.m script such that identified pre-divisional cells
%       are cropped into small tifs and labeled for single-cell comparisons
%       with Dieter's image analysis methods.

%       Labeling convention: frame number - particle ID.
%         e.g. prediv-4-2.tif refers to the particle with ID = 2 as it
%         appears in image frame number 4.



% Strategy:
%
%     0. initialize tracking and image data
%     1. isolate ellipse data from movie (stage xy) of interest
%     2. identify tracks in rejects. all others are tracked.
%     3. for each image...
%            i. load current image
%           ii. loop through predivisional cells and crop
%          iii. save cropped image
%           iv. concatenate each predivisional cell's stats
%     4. save stats in kool-kids repo!


% last edit: jen, 2020 Apr 4
% commit: crop predivisional cells and save stats


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

% 0. initialize variables for collecting stats 
stats = nan(1,4); % four columns for frame, cell ID, length and width
counter = 0;

% 1. initialize image data
conversionFactor = 6.5/60;      % scope5 Andor COSMOS = 6.5um pixels / 60x magnification
img_prefix = strcat('lb-fluc-',date,'_xy', num2str(xy), 'T'); 
img_suffix = 'XY1C1.tif';


% 2. load measured data
experimentFolder = strcat('/Users/jen/Documents/StockerLab/Data/LB/',date);
cd(experimentFolder)
filename = strcat('lb-fluc-',date,'-c123-width1p4-c4-1p7-jiggle-0p5.mat');
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
isDivision = [div_temp; 1];%length(xyData_fullOnly)];
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
names = {imgDirectory.name};
clear img_folder img_prefix img_suffix experiment newFolder img_folder


% 8. identify tracks and birth events present in each frame
totalFrames = length(imgDirectory); % total frame number
trackIDs = [];
for fr = 1:totalFrames
    
    tracksInCurrentFrame = xyData_fullOnly(xyData_fullOnly(:,16) == fr,:);   % col 16 = frame #
    trackIDs{fr,1} = tracksInCurrentFrame(:,1);                              % col 1 = TrackID, as given by tracking script ND2Proc_XY
    birthEvents{fr,1} = tracksInCurrentFrame(:,4);                           % col 4 = isDrop
    divEvents{fr,1} = isDivision(xyData_fullOnly(:,16) == fr,:);
    
end
clear fr tracksInCurrentFrame


% 9. for each image in movie, crop each newborn cell
for img = 1:length(names)
    
    % i. initialize current image
    cla
    Image = imread(names{img});
    filename = strcat('dynamicOutlines-birthEvents-xy',num2str(xy),'-frame',num2str(img),'.tif');
%     figure(1)
%     imshow(Image,'DisplayRange',[2200 4800]); %lowering right # increases num sat'd pxls
%     
    
    % ii. determine number of cells at birth in current frame
    %fr_births = birthEvents{img};
    fr_predivs = divEvents{img};
    %cells_birth = sum(fr_births);
    cells_prediv = sum(fr_predivs);
    
    
    % iii. continue to next image if no cells in either category
    if isempty(trackIDs{img}) == 1
        continue
        
    elseif cells_prediv == 0  % ii. also if no birth events 
        continue
        
    else
        
        % iii. else, isolate newborm data for each image
        is_prediv = divEvents{img,1};
        dm_currentImage = xyData_fullOnly(xyData_fullOnly(:,16) == img,:);    % col 16 = frame #
        predivData_currentImage = dm_currentImage(is_prediv == 1,:);
        
        majorAxes = getGrowthParameter(predivData_currentImage,'length');         %  lengths
        minorAxes = getGrowthParameter(predivData_currentImage,'width');          %  widths
        
        centroid_X = getGrowthParameter(predivData_currentImage,'x_pos');          % x coordinate of centroid
        centroid_Y = getGrowthParameter(predivData_currentImage,'y_pos');          % y coordinate of centroid
        angles = getGrowthParameter(predivData_currentImage,'angle');             % angle of rotation of fit ellipses

        IDs = getGrowthParameter(predivData_currentImage,'trackID');              % track ID as assigned in ND2Proc_XY
        
        
        % iv. for each cell at birth in current image, draw ellipse in Tomato
        if cells_prediv > 0
            
            for pp = 1:length(IDs)
                
                counter = counter + 1;
                currentP_b = IDs(pp); %find(IDs == birth_IDs(pp));
                
                % determine (x,y) coordinate of upper left pixel in crop
                x_coord = (centroid_X(pp)/conversionFactor) - 50;
                y_coord = (centroid_Y(pp)/conversionFactor) - 50;
                
                
                % visually confirm we indeed have the correct rectangle to crop!
                Rect = [x_coord y_coord 100 100];
%                 figure(1)
%                 hold on
%                 rectangle('Position',Rect)
                
                Image_crop = imcrop(Image,Rect);
                output = im2uint16(zeros(101,101));
                output(1:size(Image_crop,1),1:size(Image_crop,2)) = Image_crop;
                imwrite(output,strcat('prediv-', num2str(img), '-' , num2str(IDs(pp)) , '.tif'));
                
                % store stats for each prediv cell
                stats(counter,1) = img; % column 1 = frame
                stats(counter,2) = IDs(pp); % column 2 = cell ID
                stats(counter,3) = majorAxes(pp); % column 3 = length
                stats(counter,4) = minorAxes(pp); % column 4 = width
                
            end
            
        end
        clear pp color currentP_b
        
    end
    
end
cd('/Users/jen/kool-kids/')
save('check7_stats.mat','stats')


