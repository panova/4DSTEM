function [ output, Greg_stack ] = dftregistration_alignment( mask, stack_centers, STACK, usfac)

dimension = size(stack_centers);
dimension2 = size(STACK);
maskfft = fft2(mask*1);

Greg_stack = zeros(dimension2);

for i = 1:dimension(3)
    slice = stack_centers(:, :, i);
    slicefft = fft2(slice);
    
    [output] = dftregistration(maskfft, slicefft, usfac);
   
    
    diffphase = output(2);
    row_shift = output(3);
    col_shift = output(4);
    
    slice_to_shift = single(STACK(:,:, i));
    slice_to_shiftfft = fft2(slice_to_shift);
    
    if (usfac > 0),
        [nr,nc]=size(slice_to_shiftfft);
        Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
        Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
        [Nc,Nr] = meshgrid(Nc,Nr);
        Greg = slice_to_shiftfft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
        Greg = Greg*exp(1i*diffphase);
    elseif (usfac == 0)
        Greg = slice_to_shiftfft*exp(1i*diffphase);
    end
    Greg_stack(:, :, i) = abs(ifft2(Greg));
    
%     figure(2)
%     clf
%     imagesc(Greg_stack(:, :, i));
%     axis equal off
% 
end



end

