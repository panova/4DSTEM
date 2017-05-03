function [ output_args ] = Grain_stacking( stack, mycmap )
%GRAIN_STACKING Summary of this function goes here
%   Detailed explanation goes here
% 
%  figure(567); 
%  clf(); 
%  imagesc(squeeze(mean(mean(stack, 1), 2))), 
%  axis equal off;  
%  colormap(violetFire(256)) ;
%  
 
% stack = reshape(grain_array, [132 132 130*130]);

dimension = size(stack);

set_orientation = 31;
del_orientation = 10;

avg_orientations = zeros([130, 130]);

initial_grain = zeros([dimension(1), dimension(2)]);
for z = 1:dimension(3);
    for t = 1:dimension(4);
        grain = stack(:, :, z, t);
        
        vals = grain(grain>0);
        avg_orientation = mean(vals(:))
        if (avg_orientation > set_orientation - del_orientation) && (avg_orientation < set_orientation + del_orientation)
            'yes'
            initial_grain = max(initial_grain, grain);
        end
        avg_orientations(z, t) = avg_orientation;
    end
end


% figure(4567)
% clf();
% imagesc(reshape(avg_orientations, [130 130]))
% axis equal off;
% colormap(violetFire(256))

alph = (initial_grain ~= 0)*1;
figure(set_orientation)
clf();
imagesc(initial_grain, 'AlphaData', alph)
axis equal off;
colormap(mycmap)
caxis([0 180])



end

