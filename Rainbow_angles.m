function [ FFA_all, alpha ] = Rainbow_angles( peaks, realDimx, realDimy, center, mycmap )
%Peaks is an array of the found locations of the diffraction spots. It has
%the following structure:
% [ maximum number of peaks per DP 
%   2 (x and y of peak)
%   realDimx*realDimy(the total number of DPs in the STACK) ]

dimension_peaks = size(peaks);

for ID = 1:dimension_peaks(3)
    for a = 1:dimension_peaks(1)
        peak = peaks(a, 1, ID);
        if peak ==0
            peaks(a, :, ID) = NaN;
        end
    end
end

radial_distance = squeeze(sqrt( (peaks(:, 1, :) - center(1)).^2 + ...
                        (peaks(:, 2, :) - center(2)).^2  ));

b = isnan(radial_distance);
%1 pixel = 0.0410693 nm^-1
radial_distance_nm = 1./(0.0410693*radial_distance);

% figure(12);
% clf();
% imagesc(reshape(radial_distance_nm*10, [realDimx, realDimy]));
% colorbar();
% axis equal off;

size(radial_distance_nm(b==0));
mean_radial_distance = mean(radial_distance_nm(b==0));
std_radial_distance_nm = std(radial_distance_nm(b==0));

%subtracting the center to bring it to the origin:
peaks(:, 1, :) = peaks(:, 1, :) - center(1);
peaks(:, 2, :) = peaks(:, 2, :) - center(2);

intensities = peaks(:, 3, :);

FFA = rad2deg(atan(     squeeze(peaks(:, 1, :))    ./   squeeze(peaks(:, 2, :))   ));
dimension = size(FFA);
                                    
k = ~isnan(FFA);

size(k(k==1));

% figure(2);
% imagesc(-reshape(FFA, [realDimx, realDimy]));
% 
% colormap(gray);
% axis equal off;

% % For 6-fold symmetry:
% for ID = 1:dimension(2);
%     for a = 1:dimension(1);
%         ang = FFA(a, ID);
%         if ang <= -30
%             FFA(a, ID) = FFA(a, ID) + 60;
%         elseif ang >= 30;
%             FFA(a, ID) = FFA(a, ID) - 60;
%         elseif isnan(ang) == 1;
%             FFA(a, ID) = -100;
%         end
%     end
% end
% t = [0 60];

% For 180 degree symmetry: 

for ID = 1:dimension(2)
    for a = 1:dimension(1)
        ang = FFA(a, ID);
        if ang <= 0
            FFA(a, ID) = FFA(a, ID) + 180;
        elseif isnan(ang) == 1;
            FFA(a, ID) = -1;
        end
    end
end

t = [0 180];


min(FFA(:));
max(FFA(:));

% For 4-fold symmetry:
% for ID = 1:dimension(2)
%     for a = 1:dimension(1)
%         ang = FFA(a, ID);
%         if ang <= 0
%             FFA(a, ID) = FFA(a, ID) + 90;
%         elseif isnan(ang) == 1
%             FFA(a, ID) = -1;
%         end
%     end
% end
% t = [0 90];

% figure(44);
% clf();
% plot(FFA);

FFA_all = cat(3, FFA, squeeze(intensities));
% FFA_all = sort(FFA, 1);
FFA = FFA(:, :, 1);

if dimension(1) ~=1;
    'i'
    [FFAsorted, ind] = sort(FFA, 1);
    FFA = squeeze(FFAsorted(2, :))';
    dimension = size(FFA);
    k = ((FFA~=-1)*1);
%     k = ~isnan(FFA);
    size(k(k==1));

end;

%If there is only one peak: 
FFA = squeeze(reshape(FFA, [dimension(2), realDimx, realDimy]));

alpha = reshape(k, [realDimx, realDimy]);
 
%If there are more than one peak per DP; we take the brightest one:
% FFA = squeeze(reshape(FFA(1, :), [1, realDimx, realDimy]));
% alpha = reshape(k(1, :), [realDimx, realDimy]);

% alpha = alpha(realign);
jet_wrap = vertcat(jet,flipud(jet));

size(FFA_all)

figure(39);
% imagesc(squeeze(FFA));
FFA = medfilt2(squeeze(FFA));
% FFA = squeeze(FFA);
imagesc( FFA , 'AlphaData', alpha);
imagesc( medfilt2(reshape(mean(FFA_all(:, :, 1), 1), [128 128])));

colormap(mycmap);
caxis(t);
% colorbar;
axis image off





% slices = zeros(realDimx, realDimy, 13);

% ang_step = 5;
% 
% FFA_cum = zeros([realDimx, realDimy]);

% for i=-30:ang_step:30
%     i = -i;
%     FFA_mod = FFA;
% 
%     FFA_mod(FFA_mod>i+5) = 0;
%     FFA_mod(FFA_mod<i-5) = 0;
%     
%     FFa = squeeze(FFA_mod(1, :, :));
%     FFb = squeeze(FFA_mod(2, :, :));
% %     FFc = squeeze(FFA_mod(3, :, :));
%     
%     FFA_mod = FFa + FFb ;%+ FFc;
%     
%     slices(:, :, ((i+30)/ang_step)+1) = FFA_mod;
%     FFA_cum(FFA_cum==0) = FFA_mod(FFA_cum==0);
%     
% 
% end

% figure(4)
% for i= 1:2
%     hold on
%     h=imagesc(circshift(slices(:, :, i), [0, 5]));
% %     h=imagesc((slices(:, :, i)));
%     set(h,'AlphaData',0.5);
% end
% colorbar;
% colormap(gray);
% caxis([-30 30]);
% hold off
% axis equal off



end

