function [ peaks, diff_VDF ] = Iterative_DP_analysis( STACK, stack_max, center, test)
%This function takes in a STACK of DP (4D-STEM) and analyzes it slice
%by slice. 

dimension = size(STACK);
num_peaks = 1;

%For PE dataset 14:
% center = [264, 218]; 
% lattice_rad = 110;

%For PE dataset 3:
% center = [246, 209]; 
% lattice_rad = 110;

%For dataset 18
% center = [243, 220]; 

%For dataset 59
lattice_rad = 70;

%For PE dataset 14:
spot_rad = 14;

%For dataset 65:
lattice_rad = 37;
spot_rad = 7;

%For dataset 75:
lattice_rad = 70;
spot_rad = 20;

%For dataset 61:
lattice_rad = 30;
spot_rad = 7;

[column, row] = meshgrid(1:dimension(2), 1:dimension(1));

IDs = 1:dimension(3);

if test == 'True'
    IDs = [4336:4360];
%     IDs = [21190];
%     IDs = [4120];
end

diff_VDF = zeros([dimension(3), 1]);
peaks = zeros([3, dimension(3)]);

for ID = IDs;
    ID
   %Extracting the DP of interest and converting to single: 
   I = single(STACK(:, :, ID));
   I = I + abs(min(I(:)));
   I = I/stack_max;
   %Keeping original for teseting and display: 
   Ic = I;
   
   I = medfilt2(I);
   
   %masking the center spot:
   Ra = sqrt((center(1) - row).^2 + (center(2)-column).^2);
   central_masked = (lattice_rad) - Ra + 0.5;
   I(central_masked>0) = 0;
%    I(central_masked < -15-25) = 0;
   Ic(central_masked>0) = 0;
%    Ic(central_masked < -15-25) = 0;

%    figure(45)
%    clf()
%    imagesc(I);
%    axis equal off
%    caxis([0.02, 0.07])

   
   %Convolving with a pretty aggressive Gaussian:
   %For Dataset 75:
%    I = conv2(I, fspecial('gaussian', spot_rad+10, 10), 'same');
   %For dataset 61:
   I = conv2(I, fspecial('gaussian', spot_rad+10, 3), 'same');

   
%    Ic = Ic/max(Ic(:));
   
   k = sum(Ic(Ic~=0));

   diff_VDF(ID) = k;
   
   if test == 'True'
       figure(11)
       clf();
       imagesc(Ic);
       axis equal off;
       colormap(gray);

       
       figure(31)
       clf();
       imagesc((I));
       axis equal off;
       
   end
    
  
   %Finding the peak locations on the DP:    
   for i = 1:num_peaks
       Ipeak  = Peak_finder( I, 1 );
       aperture_location = [Ipeak(1, 1) Ipeak(1, 2)];
       Ra = sqrt((aperture_location(1) - row).^2 + (aperture_location(2)-column).^2);
       
       %We also mask the symmetrical peak
       Ra2 = sqrt(( (-aperture_location(1) + 2*center(1)) - row).^2 + ( (-aperture_location(2) + 2*center(2))-column).^2);
       
       masked = spot_rad - Ra + 0.5;
       masked2 = spot_rad - Ra2 + 0.5;
       
       
       spot = I(masked>0)*100; %For dataset 75 
       spot2 = I(masked2>0);
       
       if test == 'True'
           max(spot(:))
       end
              
       %For dataset 14:
       %if max(spot(:)) >= 0.0085;
       %For dataset 3:
       %if max(spot(:)) >= 0.008;
       %for dataset 18:
%        if max(spot(:)) >= 0.007;
       %for dataset 59:
%        if max(spot(:)) >= 0.038; These are all for the peak on the
%        GAUSSIAN IMAGE; for 65, we take the max of the ORIGINAL IMAGE
       %for dataset 65:
%        if max(spot(:)) > 0.05;
       %for dataset 61:
       if max(spot(:)) > 0.001;
           ff = sprintf('%3d', max(spot(:)))
           
           Ic(masked2>0) = Ic(masked2>0) + 0.05;
           Ic(masked>0) = Ic(masked>0) + 0.03;

           I(masked>0) = 0;
           I(masked2>0) = 0;
           
           peaks(:, ID) = Ipeak;
           
       end

        
   end
   
   if test == 'True'
      
       figure(121)
       clf();
       imagesc(Ic);
       axis equal off;
       colormap(gray);
       title(ff);

   end
   
   %Erasing the DP completely, just to make sure that there isn't any 
   %cummulating arrays that would overwhelm the memory:
   I = []; 
   
end

% diff_VDF = reshape(diff_VDF, [129, 128]);

% figure(7)
% clf()
% imagesc(diff_VDF);
% axis equal off;



end

