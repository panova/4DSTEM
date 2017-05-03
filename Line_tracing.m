function [ All_coords, Intensities ] = Line_tracing( FFA, colormap )
%LINE_TRACING Summary of this function goes here
%   Detailed explanation goes here

dimension = size(FFA);
FFA = reshape(FFA, [dimension(1), 128, 128, 2]);
dimension = size(FFA);
%Dimensions of FFA: num peaks, real dimensions x, then y, 
%then THETA (stored in :, :, :, 1) and INTENSITIES (stored in :, :, :, 2)
FFA(:, :, :, 2) = FFA(:, :, :, 2)./max(FFA(:, :, :, 2));

%Row coordinates first, columns second; num peaks; real space dims
All_coords = zeros([2, 2, 4, 128, 128]);
Intensities = zeros([4, 128, 128]);

for globalC = 1:128; %55
    globalC
    %     for R = 2:realDimX-1
    %This is the plotted Y on imagesc
    for globalR = 1:128 %82
%         pause(0.05);
    
        globalR
        peaks = FFA(:, globalR, globalC, :);
        
        for peakID = 1:4
            peak = squeeze(peaks(peakID, 1));
            intensity = squeeze(peaks(peakID, 2));
            if peak == -1;
                'There is no diffraction here!' ;              
                continue
            end     
            
            m = tan(degtorad(peak));
            b = -globalR - m.*globalC;
            x = [globalC - 0.5, globalC+0.5];
            y = m.*x + b;
            
            if min(y(:)) <= -globalR-0.5;
                y = [-globalR+0.5, -globalR-0.5];
                x = (y-b)./m;
            end
            
            All_coords(1, :, peakID, globalR, globalC) = x;
            All_coords(2, :, peakID, globalR, globalC) = y; 
            Intensities(peakID, globalR, globalC) = intensity;
            
%             figure(34);
%             hold on;
%             plot(x, y, 'k', 'Color', [0 0 0] + intensity);
%             axis equal off;
%             hold off;
%             drawnow;
        end
%         drawnow;
    
    end

end




% Intensities = zeros([3, 128*128*4]) + squeeze(reshape(Intensities, [1, 128*128*4]));

% figure(3546);
% clf();
% hold on;
% plot(reshape(All_coords(1, :, :, :, :), [2, 128*128*4]), reshape(All_coords(2, :, :, :, :), [2, 128*128*4]), 'k');
% axis equal off; 
% hold off;

All_coordsx = reshape(All_coords(1, :, :, :, :), [2, 128*128*4]);
All_coordsy = reshape(All_coords(2, :, :, :, :), [2, 128*128*4]);
base_alpha = zeros([4, 128*128*4]);
base_alpha(4, :) = squeeze(reshape(Intensities, [1, 128*128*4]));

%rescaling the colormap to match the thetas: 
L = size(colormap, 1);
thetas = FFA(:, :, :, 1);
thetas = squeeze(reshape(thetas, [1, 128*128*4]));

upcolormap1 = interp1(1:(182/64):182, colormap(:, 1), 1:180);
upcolormap2 = interp1(1:(182/64):182, colormap(:, 2), 1:180);
upcolormap3 = interp1(1:(182/64):182, colormap(:, 3), 1:180);

size(upcolormap1)
% %For some reason the interpolation does not go all the way to the end if it stops at 180,
% %so it appends some "NaN"s at the end. Those are erased. 
% upcolormap1(end-1:end) = upcolormap1(end-2);

thetas(thetas < 1) = 1;

base_alpha(1, :) = squeeze(upcolormap1(round(thetas(:))));
base_alpha(2, :) = squeeze(upcolormap2(round(thetas(:))));
base_alpha(3, :) = squeeze(upcolormap3(round(thetas(:)))); 

figure(3546);
clf();
hold on;
for t = 1:128*128*4
    t
% plot(All_coordsx(:, t), All_coordsy(:, t), 'Color', Intensities(:, t));
plot(All_coordsx(:, t), All_coordsy(:, t), 'Color', base_alpha(:, t));

end
axis equal off; 
hold off;

print -painters -dpdf -r600 TCDIO4_Linear_raster_map_FULLALPHACOLOR.pdf

end


