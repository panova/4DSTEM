function [ Greg_stack ] = dft_registration_SHIFTING( STACK, output )
%DFT_REGISTRATION_SHIFTING Summary of this function goes here
%   Detailed explanation goes here

dimension = size(STACK);

format = 'uint16=>uint16';
Greg_stack = zeros(dimension, 'uint16');

for i = 1:dimension(3)
        slice = STACK(:, :, i);
        slicefft = fft2(slice);
        
        diffphase = output(i, 2);
        row_shift = output(i, 3);
        col_shift = output(i, 4);
        
        [nr,nc]=size(slicefft);
        Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
        Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
        [Nc,Nr] = meshgrid(Nc,Nr);
        Greg = slicefft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
        Greg = Greg*exp(1i*diffphase);
        
        Greg_stack(:, :, i) = abs(ifft2(Greg));
        
        %     figure(2)
        %     clf
        %     imagesc(Greg_stack(:, :, i, i2));
        %     axis equal off

end


end

