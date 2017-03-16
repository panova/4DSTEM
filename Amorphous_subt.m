function [ subSTACK,  amorphous_map ] = Amorphous_subt( STACK, center, average, lattice_rad, realDimx, realDimy )
%This function takes a slice from a stack identified by its ID number within the stack and subtracts the average amorphous halo from it. 

dimension = size(STACK);
temp_radius = lattice_rad;
% center = [dimension(1)/2, dimension(2)/2];
% center = [266+8, 148+91];

%creation of the amorphous signal aperture (taken to be close to the
%beamstop)
[Ya, Xa] = meshgrid(1:dimension(2), 1:dimension(1));
Ra = sqrt ( (center(1) - Xa).^2     +      (center(2) - Ya).^2 )    ;
aperture = temp_radius - Ra;
aperture(aperture<0) = 0;
aperture(aperture>0) = 1;
aperture = (aperture~=0)*1;

figure(27);
clf();
imagesc(aperture.*STACK(:, :, 1808));
axis equal off;
drawnow


%This is the actual virtual dark field generated from the "amorphous"
%signal located within the aperture created above. 
map = VIRTUAL_DF(STACK, aperture, realDimx, realDimy);
amorphous_map = map/max(map(:));

figure(23);
clf();
imagesc(amorphous_map);
axis equal off;
drawnow;

%Each real space pixel of the amorphous VDF is replicated to be the size of
%the DP it is supposed to be subtracted from. 
map = reshape(map, [1, 1, dimension(3)]);

size(map(:))

map = single(map/max(map(:)));
map = repmat(map, [dimension(1), dimension(2), 1]);
size(map(:))

averages = repmat(average, [1, 1, dimension(3)]).*map*3;
subSTACK = STACK - averages;

subSTACK(subSTACK <0) = 0;

subSTACK = 255*(subSTACK/max(subSTACK(:)));
end

