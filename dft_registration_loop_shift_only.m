function [ output ] = dft_registration_loop_shift_only( stack_centers, usfac )
%DFT_REGISTRATION_LOOP_SHIFT_ONLY Summary of this function goes here
%   Detailed explanation goes here

dimension = size(stack_centers);
mask = stack_centers(:, :, 1);
maskfft = fft2(mask*1);

output = zeros([dimension(3), 4]);

for i = 1:dimension(3)
    slice = stack_centers(:, :, i);
    slicefft = fft2(slice);
    
    output(i, :) = dftregistration(maskfft, slicefft, usfac);

end

end

