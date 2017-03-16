function [ filmarray ] = Film_import( )
%MATRIX_FROM_PNGS Summary of this function goes here
%   Detailed explanation goes here


% myFolder = '\\tsclient\COLIN PC REDIRECT\Winey PE 2016\20160321 Winey HDPE xtal#1 Movie data\28Diff during STEM acquis dwell=10us Cl=380 spot11 pixel size=4p7nm 128x128_28';
% myFolder = '\\tsclient\COLIN PC REDIRECT\Winey PE 2016\20160421 For OUliana STEM movies 29 to 41\32';
myFolder = 'C:\RESEAARSH\Time Constant Movies PE\40';


if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

filePattern = fullfile(myFolder, '*.png');
jpegFiles = dir(filePattern);
length(jpegFiles)

filmarray = [];
for k = 1:length(jpegFiles)
  baseFileName = jpegFiles(k).name;
  
  fullFileName = fullfile(myFolder, baseFileName);
  fprintf(1, 'Now reading %s\n', fullFileName);
  imageArray = rgb2gray(imread(fullFileName));
  filmarray = cat(3, filmarray, imageArray);
  %imshow(imageArray);  % Display image.
  %drawnow; % Force display to update immediately.
end

end

