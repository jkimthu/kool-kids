%% check five: cropped single cells for Dieter


% Goal: Modify check1.m script such that identified newborns are cropped
%       into small tifs and labeled for single-cell comparisons with
%       Dieter's image analysis methods.

%       Labeling convention: frame number - particle ID.
%         e.g. crop-4-2.tif refers to the particle with ID = 2 as it
%         appears in image frame number 4.



% Strategy:
%
%     0. initialize tracking and image data
%     1. isolate ellipse data from movie (stage xy) of interest
%     2. identify tracks in rejects. all others are tracked.
%     3. for each image...
%            i. load current image
%           ii. loop through newborn cells and crop
%          iii. save cropped image
%     4. woohoo!


% last edit: jen, 2019 Mar 29
% commit: cropped newborns for dieter


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
names = {imgDirectory.name};
clear img_folder img_prefix img_suffix experiment newFolder img_folder



% 8. identify tracks and birth events present in each frame
totalFrames = length(imgDirectory); % total frame number
trackIDs = [];
for fr = 1:totalFrames
    
    tracksInCurrentFrame = xyData_fullOnly(xyData_fullOnly(:,16) == fr,:);   % col 16 = frame #
    trackIDs{fr,1} = tracksInCurrentFrame(:,1);                              % col 1 = TrackID, as given by tracking script ND2Proc_XY
    birthEvents{fr,1} = tracksInCurrentFrame(:,4);                           % col 4 = isDrop
    
end
clear fr tracksInCurrentFrame


% 9. for each image in movie, crop each newborn cell
for img = 1:length(names)
    
    % i. initialize current image
    cla
    Image = imread(names{img});
    filename = strcat('dynamicOutlines-birthEvents-xy',num2str(xy),'-frame',num2str(img),'.tif');
    %figure(1)
    %imshow(Image,'DisplayRange',[2200 4800]); %lowering right # increases num sat'd pxls
    
    
    % ii. determine number of cells at birth in current frame
    fr_births = birthEvents{img};
    cells_birth = sum(fr_births);
    
    
    % iii. continue to next image if no cells in either category
    if isempty(trackIDs{img}) == 1
        continue
        
    elseif cells_birth == 0  % ii. also if no birth events 
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
                imwrite(output,strcat('crop-', num2str(img), '-' , num2str(birth_IDs(pp)) , '.tif'));
                
            end
            
        end
        clear pp color birth_IDs currentP_b
        
    end
    
end



