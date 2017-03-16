function [ VDF_STACK ] = VirtualDF_stack_from_rotating_aperture( STACK, mask_STACK )
%BE AWARE: dimensions of the template stack are (angleID, DPx, DPy) - which
%is NOT THE SAME as the stack, which is (DPx, DPy, sliceID). 
%mask_STACK is the stack of templates with the rotating aperture. 

dimension = size(mask_STACK);
% realDimx = 128;
% realDimy = 128;

realDimx = 50;
realDimy = 50;

VDF_STACK = zeros([realDimx, realDimy, dimension(1)]);
size(VDF_STACK)

for i = 1:dimension(1);
    VDF_STACK(:, :, i) = VIRTUAL_DF( STACK, squeeze(mask_STACK(i, :, :)), realDimx, realDimy );

end


end

