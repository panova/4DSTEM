function [ subSTACK ] = Rainbow_prep( STACK )
%RAINBOW_PREP Summary of this function goes here
%   Detailed explanation goes here

dimension = size(STACK);

%Dimensions of the capture window in real space (field of view, in pixels)
realDimx = 128;
realDimy = 128;

[ XX, YY ] = Linear_center_shifts( realDimx, realDimy );

STACK = shifting_stack_linear_XXYY( STACK, XX, YY );

% STACK = permute(STACK, [2, 1, 3]);

%Converting to single:
STACK = single(STACK);

mean_STACK = mean(STACK, 3);

center = [dimension(1)/2 dimension(2)/2];

% center = [125, 125]; %for A1, A2
%The radius in pixels of the diffraction ring to analyze
lattice_rad = 65;

%Subtracting the amorphous background
[ subSTACK,  ~ ] = Amorphous_subt( STACK, mean_STACK, lattice_rad, realDimx, realDimy );

%After the amorphous subtraction, we threshold all the negative values to
%zero. 
%NOTE: The amorphous background subtraction is set up to be based on
%"amorphous" signal close to the central beam, which means that the
%subtraction may be too harsh at times. Adjustment of the multiplication
%constant of the average to be subtracted may have to be adjusted. 
% subSTACK(subSTACK<=0) = 0;

mean_subSTACK = mean(subSTACK, 3);

figure(15)
imagesc(mean_subSTACK);
axis equal off




end

