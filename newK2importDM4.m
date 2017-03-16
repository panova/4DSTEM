function [stack3D,xv,yv] = newK2importDM4(topDir,...
    minuteRange,secondRange,frameRange)
% Import "DM4" data from K2 camera into 3D dataset (x,y,frame)
%  - topDir = the top director that contains the Hour_00 directory just below it
%  - secondRange = second directories to read in (0:1 for example to read second_00 and second_01 dirs)
%  - frameRange = Optional frames to read in each second directory. Start at 0, max is 399. Leave blank to read all frames in each second directory
%  - Based on Colin's quickImport02.m.
%  - Added "minuteRange" -Colin

cropSize = [512 512];
cropSize = [256 256];
% cropSize = [573 573]; %For dataset 75

% cropSize = [1920 1792/2-1];
cropOrigin = [1 1];

% cropSize = [700 760];
% cropOrigin = [631 517];


% Data structure information
dirBase = fullfile(topDir,'Hour_00');
dd = dir(fullfile(dirBase,'Minute_00','Second_00','*.dm4'));
fileBase = strsplit(dd(1).name,'Hour_');
fileBase = char(fileBase(1));
offset = 0; % in numerical indices!!
%format = 'uint16=>unit16'; %these are determined later
%xySize = [1920 1792]; %these are determined later
fileStem = '.dm4';

dirHour = 'Hour_00_';
% dirMin = 'Minute_00_';

if nargin == 3
    frameRange = 0:length(dd)-1;
    disp('Reading all files in each seconds directory')
end

%Test if directory contains at least as many files as requested
if length(frameRange) > length(dd)
    error('requested more files than available in the directoy');
end

%Open the first image to get the file properties
dirMin0 = ['Minute_', sprintf('%02d',minuteRange(1))];
dirSec0 = ['Second_', sprintf('%02d',secondRange(1))];
frameNum0 = ['Frame_', sprintf('%04d',frameRange(1))];
firstFileName = [fileBase,dirHour,dirMin0,'_',dirSec0,'_',frameNum0,fileStem];
fname = fullfile(dirBase,dirMin0,dirSec0,firstFileName);
im0 = dm4Reader(fname);
try
    format = [class(im0.data),'=>',class(im0.data)];
    fc = class(im0.data);
    xySize = size(im0.data);
end
try
    format = [class(im0.image),'=>',class(im0.image)];
    xySize = size(im0.image);
    fc = class(im0.image);
end

%Cropping variables
% cropSize = xySize;
% cropOrigin = [1 1];
% cropSize = [xySize(1) 142];
% cropOrigin = [1 187];


% Main loop
inds = offset+(1:(prod(xySize)));
xv = cropOrigin(1)-1+(1:cropSize(1));
yv = cropOrigin(2)-1+(1:cropSize(2));
% xv = cropOrigin(1)-1+((-cropSize(1)/2):(cropSize(1)/2-1));
% yv = cropOrigin(2)-1+((-cropSize(2)/2):(cropSize(2)/2-1));
Nstack = [cropSize(1:2) length(minuteRange)*length(secondRange)*length(frameRange)];
% stack3D = zeros(Nstack(1),Nstack(2),Nstack(3),class(im0.data));
stack3D = zeros(Nstack(1),Nstack(2),Nstack(3),fc);


numFrames = length(frameRange);
numSeconds = length(secondRange);
numMinutes = length(minuteRange);

% progressbar(0,2); %initialize progressbar
errorCount = 0;

%Time the import
% tic
count = 0;  % number of frames!
for a0 = 1:numMinutes
    dirMin = ['Minute_', sprintf('%02d',minuteRange(a0))];
    
    for a1 = 1:numSeconds
        dirSec = ['Second_', sprintf('%02d',secondRange(a1))]
        
        %Build filename and load file
        for a2 = 1:numFrames
            
            % Build file name and directory name
            frameNum = ['Frame_', sprintf('%04d',frameRange(a2))];
            fileBuild = [fileBase,dirHour,dirMin,'_',dirSec,'_',frameNum,fileStem];
            
            % Open the DM4 file
            fname = fullfile(dirBase,dirMin,dirSec,fileBuild);
            if exist(fname,'file')
                
                fid = fopen(fname,'rb');
                fseek(fid,im0.binaryOffset(2),-1); %skip to the start of the data
                
                if(fid ~= -1) %check if file exists
                    dataStream = fread(fid,xySize(1)*xySize(2),format);
                    fclose(fid);
                    
                    % Reshape image data and crop ROI
                    frame = reshape(dataStream(inds),xySize);
                    count = count + 1;
                    stack3D(:,:,count) = ...
                        frame(xv,yv);
%                     stack3D(:,:,count) = ...
%                         frame(:,:)';
                    
                else %if file does not exist
                    errorCount = errorCount + 1;
                    fname1 = strsplit(fname,'\');
                    notExist{errorCount} = fname1(end);
                    if errorCount == 1
                        disp('a file not found. Wait for error.')
                    end
                end
            end
        end
        
        % Progressbar
        comp = (a1 / numSeconds ...
            + a0 - 1) /numMinutes;
        progressbar(comp,2);
    end
end
% toc
% progressbar(1,2); %end progressbar if it still is open

% % Remove excess stack frames
% if size(stack3D,3) > count
%     stack3D(:,:,(count+1):end) = [];
% end

if(errorCount > 0)
    errordlg('See command window for non existent file names.')
    disp(notExist')
end

end