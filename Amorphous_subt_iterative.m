function [ STACK,  amorphous_map ] = Amorphous_subt_iterative( STACK, center, average, lattice_rad, realDimx, realDimy )
%This function takes a slice from a stack identified by its ID number 
% within the stack and subtracts the average amorphous halo from it. 

dimension = size(STACK);
temp_radius = lattice_rad;
% center = [dimension(1)/2, dimension(2)/2];


%creation of the amorphous signal aperture (taken to be close to the
%beamstop)
[Ya, Xa] = meshgrid(1:dimension(2), 1:dimension(1));
Ra = sqrt ( (center(1) - Xa).^2     +      (center(2) - Ya).^2 )    ;
apertureBeamstop = 50 - Ra;
apertureBeamstop(apertureBeamstop<0) = 0;
apertureBeamstop(apertureBeamstop>0) = 1;
apertureBeamstop = (apertureBeamstop==0)*1;


%creation of the amorphous signal aperture (taken to be close to the
%beamstop)
[Ya, Xa] = meshgrid(1:dimension(2), 1:dimension(1));
Ra = sqrt ( (center(1) - Xa).^2     +      (center(2) - Ya).^2 )    ;
aperture = temp_radius - Ra;
aperture(aperture<0) = 0;
aperture(aperture>0) = 1;
aperture = (aperture==0)*1;

figure(234);
clf();
imagesc(STACK(:, :, 12821).*apertureBeamstop);
axis equal off;
drawnow;

figure(45);
clf();
imagesc(STACK(:, :, 12821).*aperture);
axis equal off;
drawnow;

for i=1:dimension(3)
    i
    STACK(:, :, i) = STACK(:, :, i).*apertureBeamstop;
end

%This is the actual virtual dark field generated from the "amorphous"
%signal located within the aperture created above. 
map = VDF_iterative(STACK, aperture);
map = reshape(map, [realDimx, realDimy]);
amorphous_map = map/max(map(:));

%Each real space pixel of the amorphous VDF is replicated to be the size of
%the DP it is supposed to be subtracted from. 
map = reshape(map, [1, 1, dimension(3)]);
map = single(map/max(map(:)));

for i=1:dimension(3)
    i
    slice = single(STACK(:, :, i));
    m = squeeze(map(:, :, i));
    slice = slice - average*m*2;
    STACK(:, :, i) = slice;
    
end



figure(34);
clf();
imagesc(amorphous_map);
axis equal off;
drawnow;

STACK(STACK<0) = 0;
STACK = 255*(STACK/max(STACK(:)));

end

