function [ Peaks_all ] = Rainbow_map( subSTACK, center, lattice_rad, Test )
%This function takes in an angle step (10 for every 10 degrees, etc) and
%makes a rainbow map where the real-space image is colored according to the
%angle of orientation of the crystallites. This is based on the rotating
%virtual dark field technique. 
tic
dimension = size(subSTACK);

subSTACK = subSTACK./max(subSTACK(:));

% subSTACK = subSTACK(1:round(dimension(1)/2), :, :);
% aperture = aperture(1:round(dimension(1)/2), :, :);
% dilbeam = dilbeam(1:round(dimension(1)/2), :, :);
% dimension = size(subSTACK);

%Dimensions of the capture window in real space (field of view, in pixels)
realDimx = 128;
realDimy = 128;

% [ XX, YY ] = Linear_center_shifts( realDimx, realDimy );

% STACK = shifting_stack_linear_XXYY( STACK, XX, YY );

% STACK = permute(STACK, [2, 1, 3]);

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

%Subtracting the amorphous background
% [ subSTACK,  ~ ] = Amorphous_subt( STACK, mean_STACK, lattice_rad, realDimx, realDimy );

%After the amorphous subtraction, we threshold all the negative values to
%zero. 
%NOTE: The amorphous background subtraction is set up to be based on
%"amorphous" signal close to the central beam, which means that the
%subtraction may be too harsh at times. Adjustment of the multiplication
%constant of the average to be subtracted may have to be adjusted. 
% subSTACK(subSTACK<=0) = 0;
% 
% mean_subSTACK = mean(subSTACK, 3);

% figure(5)
% imagesc(mean_subSTACK);
% axis equal off

num_peaks = 1;

Peaks_all = zeros([num_peaks 2 dimension(3)]);

[column, row] = meshgrid(1:dimension(2), 1:dimension(1));

aesthetic_aperture_rad = 12; % Equal to 6 for D2; This is for display purposes and testing. This is the radius of the expected diffraction spot. 
masking_aperture_rad = aesthetic_aperture_rad*3; %When looking for peaks on the DP, the peaks found will be masked with an aperture of this radius

% middle_masked = zeros([dimension(1), dimension(2), dimension(3)]);

Ra = sqrt((center(1) - row).^2 + (center(2)-column).^2);
central_masked = (lattice_rad/1.3) - Ra + 0.5;
aperture = central_masked;
aperture(aperture>0) = 0;
aperture(aperture<0) = 1;
aperture = (aperture~=0)*1;

figure(145);
clf();
imagesc(aperture);
axis equal off;

% for t = 1:dimension(3)
%     t
%     I = aperture.*subSTACK(:, :, t);
%     subSTACK(:, :, t) = I;
% end;

subSTACK = subSTACK./max(subSTACK(:));

%We now step over each DP:
% for ID = 1620:dimension(3);
% for ID = 48:51%dimension(3);
for ID = 1:dimension(3);
    ID
    I = subSTACK(:, :, ID);
    I = single(I);
    Ik = I;

    %Deleting all pixels below 10% of the maximum brightness
    I(I<0.02) = 0;

    I = medfilt2(I);

    Im = I;


%     Central aperture masking the bright central beam:
    I = aperture.*I;
    I(central_masked>0) = 0;
    Im(central_masked>0) = 0;
%     I(central_masked<-aesthetic_aperture_rad-25) = 0;
    
%     I = I.^2;
%     I = I.*aperture.*dilbeam;
%     I = I.*aperture;
    
%     Im(central_masked<-aesthetic_aperture_rad-25) = 0;
%     middle_masked(:, :, ID) = I;
    I = Gaussian_blurr(I, 15, 1);
%     Im = I;

    if Test == 1
        figure(45);
        clf();
        imagesc(Ik);
        title(ID)
        axis equal off;
%         caxis([0 0.007]);
        drawnow();
    end
    
    Ipeaks = zeros(num_peaks, 2);
        
    for i = 1:num_peaks;
        Ipeak  = Peak_finder( I, 1 );
        aperture_location = [Ipeak(1, 1) Ipeak(1, 2)];
        Ra = sqrt((aperture_location(1) - row).^2 + (aperture_location(2)-column).^2);
        %We also mask the symmetrical peak
        Ra2 = sqrt(( (-aperture_location(1) + 2*center(1)) - row).^2 + ( (-aperture_location(2) + 2*center(2))-column).^2);

        masked = masking_aperture_rad - Ra + 0.5;
        masked2 = masking_aperture_rad - Ra2 + 0.5;
        
        aesthetic = aesthetic_aperture_rad - Ra + 0.5;
        aesthetic2 = aesthetic_aperture_rad - Ra2 + 0.5;
        
%         I = I/max(I(:));
        peak_signal = sum(I(aesthetic>0)) + sum(I(aesthetic2>0));
        st = std(I(aesthetic>0));
        
%         if peak_signal>1e4
%         if peak_signal>-700 % for E1
        if peak_signal>0
            
%             if st>20    
            if st>0.             

                peak_signal = sum(I(aesthetic>0)) + sum(I(aesthetic2>0));
                
                I(masked>0) = 0;
                I(masked2>0) = 0;
                
                Im(aesthetic>0) =  Im(aesthetic>0)+ 0.1;
                Im(aesthetic2>0) =  Im(aesthetic2>0)+ 0.05;
                
                Ipeaks(i, :) = Ipeak(1:2);
            else
                peak_signal = 0;
                Ipeaks(i, :) = [0, 0];
            end
        else
            peak_signal = 0;
            Ipeaks(i, :) = [0, 0];
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

