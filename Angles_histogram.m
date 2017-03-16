function [ FFA ] = Angles_histogram( Peaks_allPE )
%ANGLES_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here

FA = squeeze(Peaks_allPE(1, :, :));
% rFA = FA(1, :) - 227.1; %For PE
% cFA = FA(2, :) - 209.1; %For PE

rFA = FA(1, :) - 274;
cFA = FA(2, :) - 260;

FFA = rad2deg(atan(cFA./rFA));
dimension = size(FFA(:));

for ID = 1:dimension(1);
   ang = FFA(ID);
   if ang <= -30
       FFA(ID) = FFA(ID) + 60;
   elseif ang >= 30;
       FFA(ID) = FFA(ID) - 60;
   end
end

% FFA = reshape(FFA, [129 128]);
FFA = reshape(FFA, [50 50]);

min(FFA(:))
max(FFA(:))

figure(23);
clf();
imagesc(FFA);
colormap(jet)
colorbar();
axis equal off;

Peaks_allPE = reshape(Peaks_allPE, [], 2);

RPeaks_allPE = (Peaks_allPE(:, 1));
CPeaks_allPE = (Peaks_allPE(:, 2));

Mask = (RPeaks_allPE+CPeaks_allPE ~= 0);

RPeaks_allPE = RPeaks_allPE(Mask ~= 0) - 227.1;
CPeaks_allPE = CPeaks_allPE(Mask ~= 0) - 209.1;

angles = rad2deg(atan(CPeaks_allPE./RPeaks_allPE));


figure(2); clf(); 
hist(angles, 100);
% h.Normalization = 'countdensity';




end

