function [ Peaks_all ] = Rainbow_mapPC( subSTACK, center, lattice_rad, Test )
%This function takes in an angle step (10 for every 10 degrees, etc) and
%makes a rainbow map where the real-space image is colored according to the
%angle of orientation of the crystallites. This is based on the rotating
%virtual dark field technique. 
tic
dimension = size(subSTACK);

% subSTACK = subSTACK./max(subSTACK(:));

% subSTACK = subSTACK(1:round(dimension(1)/2), :, :);
% aperture = aperture(1:round(dimension(1)/2), :, :);
% dilbeam = dilbeam(1:round(dimension(1)/2), :, :);

%Dimensions of the capture window in real space (field of view, in pixels)
realDimx = 128;
realDimy = 128;

%Converting to single:
% STACK = single(STACK);

% mean_STACK = mean(STACK, 3);

% center = [123, 130]; % for D2
% center = [113, 114]; % for C1
% center = [131, 155]; % for E1
% center = [131, 141]; % for E1 halved
% center = [266+8, 148+91];


% center = [125, 125]; %for A1, A2
%The radius in pixels of the diffraction ring to analyze for D2
% lattice_rad = 64;
% lattice_rad = 60;

%After the amorphous subtraction, we threshold all the negative values to
%zero. 
%NOTE: The amorphous background subtraction is set up to be based on
%"amorphous" signal close to the central beam, which means that the
%subtraction may be too harsh at times. Adjustment of the multiplication
%constant of the average to be subtracted may have to be adjusted. 
% subSTACK(subSTACK<=0) = 0;
% This is all performed by the Amorphous_subt_iterative.m code


% mean_subSTACK = mean(subSTACK, 3);

% figure(5)
% imagesc(mean_subSTACK);
% axis equal off

num_peaks = 4;

Peaks_all = zeros([num_peaks 3 dimension(3)]);

[column, row] = meshgrid(1:dimension(2), 1:dimension(1));



aesthetic_aperture_rad = 12; % Equal to 6 for D2; This is for display purposes and testing. This is the radius of the expected diffraction spot. 
% masking_aperture_rad = aesthetic_aperture_rad*3; %When looking for peaks on the DP, the peaks found will be masked with an aperture of this radius
masking_aperture_rad = aesthetic_aperture_rad*2; %When looking for peaks on the DP, the peaks found will be masked with an aperture of this radius

% middle_masked = zeros([dimension(1), dimension(2), dimension(3)]);

%Creating the aperture mask that masks away the remnants of the central beam: (may not be subsequently used for large
%datasets)
Ra = sqrt((center(1) - row).^2 + (center(2)-column).^2);
central_masked = (lattice_rad/1.3) - Ra + 0.5;
aperture = central_masked;
aperture(aperture>0) = 0;
aperture(aperture<0) = 1;
aperture = (aperture~=0)*1;
aperture = Gaussian_blurr(aperture, 20, 5);

%Checking the location of the aperture:
% figure(145);
% clf();
% imagesc(aperture);
% axis equal off;

subSTACK = subSTACK./max(subSTACK(:));

%We now step over each DP:
% for ID = 1620:dimension(3);
% for ID = 152%dimension(3);
for ID = 1:dimension(3);
% for ID = 16384;
    ID
    I = subSTACK(:, :, ID);
    I = single(I);

    if Test == 1
        figure(495);
        clf();
        imagesc(I);
        title(ID)
        axis equal off;
        drawnow();
        
    end;
        
    %Deleting all pixels below 10% (or 20%, depending on dataset) of the maximum brightness
    I(I>0.01) = 0.01;

    I = medfilt2(I);
    Im = I;

%     Central aperture masking the bright central beam:
    I = aperture.*I;
%     I(central_masked>0) = 0;
%     Im(central_masked>0) = 0;
%     I(central_masked<-aesthetic_aperture_rad-25) = 0;
%     Im(central_masked<-aesthetic_aperture_rad-25) = 0;
%     middle_masked(:, :, ID) = I;

    if Test == 1
        figure(8);
        clf();
        imagesc(I);
        title(ID)
        axis equal off;
        %         caxis([0 0.007]);
        drawnow();
        
    end;

    PC = Phase_corellation(I, aesthetic_aperture_rad);
%     PC(central_masked>0) = 0;
    PC = PC.*aperture;
    PC = Gaussian_blurr(PC, aesthetic_aperture_rad, 1);
%     PC = medfilt2(PC);

%     Im = I;

    if Test == 1
        figure(45);
        clf();
        imagesc(I);
        title(ID)
        axis equal off;
%         caxis([0 0.007]);
        drawnow();
        
        figure(47);
        clf();
        imagesc(PC);
        title('Phase Correlation')
        axis equal off;
%         caxis([0 0.007]);
        drawnow();
    end
    
    
    
    Ipeaks = zeros(num_peaks, 3);
        
    for i = 1:num_peaks;
        Ipeak  = Peak_finder( PC, 1 );
        aperture_location = [Ipeak(1, 1) Ipeak(1, 2)];
        Ra = sqrt((aperture_location(1) - row).^2 + (aperture_location(2)-column).^2);
        %We also mask the symmetrical peak
        Ra2 = sqrt(( (-aperture_location(1) + 2*center(1)) - row).^2 + ( (-aperture_location(2) + 2*center(2))-column).^2);

        masked = masking_aperture_rad - Ra + 0.5;
        masked2 = masking_aperture_rad - Ra2 + 0.5;
        
        aesthetic = aesthetic_aperture_rad - Ra + 0.5;
        aesthetic2 = aesthetic_aperture_rad - Ra2 + 0.5;
        
%         peak_signal = sum(I(aesthetic>0)) + sum(I(aesthetic2>0));
        peak_signal = sum(I(aesthetic>0));

        st = std(I(aesthetic>0));
        
%         if peak_signal>1e4
%         if peak_signal>-700 % for E1
        if peak_signal>1
            
%             if st>20    
            if st>0.             

                peak_signal = sum(I(aesthetic>0));
                
                I(masked>0) = 0;
                PC(masked>0) = 0;

%                 I(masked2>0) = 0;
%                 PC(masked2>0) = 0;

                Im(aesthetic>0) =  Im(aesthetic>0)+ 0.04;
%                 Im = Im.*aesthetic;
%                 Im(aesthetic2>0) =  Im(aesthetic2>0)+ 0.05;
                
                Ipeaks(i, :) = Ipeak(1:3);
                Ipeaks(i, 3) = peak_signal;
            else
                peak_signal = 0;
                Ipeaks(i, :) = [0, 0, 0];
            end
        else
            peak_signal = 0;
            Ipeaks(i, :) = [0, 0, 0];
        end
        
    end
    Peaks_all(:, :, ID) = Ipeaks;
    
%     figure(1)
%     imagesc(Im.*aperture.*dilbeam);
%     axis equal off;
%     

if Test == 1
    figure(455);
    clf();
    imagesc(Im);
    axis equal off;
end
%     

% pause(0.2);
  
end



toc

end

