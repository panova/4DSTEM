function [ VDF ] = VDF_iterative( STACK, aperture )
%VDF_ITERATIVE Summary of this function goes here
%   Detailed explanation goes here


dimension = size(STACK); 

VDF = zeros([dimension(3), 1]);

for i = 1:dimension(3);
    i
   slice = single(STACK(:, :, i));
   slice = slice.*aperture;
   k = sum(slice(:));
   VDF(i) = k;
end



end

