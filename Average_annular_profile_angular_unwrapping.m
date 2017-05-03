function [ sum_I ] = Average_annular_profile_angular_unwrapping( STACK, realDimX, realDimY, center, mycmap )
%AVERAGE_ANNULAR_PROFILE_ANGULAR_UNWRAPPING Summary of this function goes here
%   Detailed explanation goes here


stack_dimension = size(STACK);
% padd = 3;

lattice_rad = 76;
spot_rad = 12;

% mean_stack = mean(STACK, 3);

STACK = reshape(STACK, [stack_dimension(1), stack_dimension(2), realDimX, realDimY]);
% STACK = padarray(STACK, [0, 0, padd, padd], 'both');

stack_dimension = size(STACK);

[Ya, Xa] = meshgrid(1:stack_dimension(2), 1:stack_dimension(1));
Ra = sqrt((center(1) - Xa).^2 + (center(2)-Ya).^2);
mask = (lattice_rad) - Ra;
mask2 = (lattice_rad + 3*spot_rad) - Ra;
mask(mask>=0) = 1;
mask(mask<0) = 0;
mask2(mask2>=0) = 1;
mask2(mask2<0) = 0;
aperture = abs((mask - mask2));

sum_I = zeros([5000, 1]);

for R = 2:(stack_dimension(3)-1)
    
    R
    for C = 2:(stack_dimension(4)-1)    
                DP = STACK(:, :, R, C);       
        
        DP = DP.*aperture;
        
        [rows, cols, I] = find(DP);
        rows = rows - center(1);
        cols = cols - center(2);
%         theta = rad2deg(atan(     rows    ./   cols   ));
        theta = rad2deg(atan2(rows, cols));
        pts = sortrows(cat(2, theta, I));
        
        theta = pts(:, 1);
        I = pts(:, 2);

        [theta, ia, ~] = unique(theta, 'stable');
        I = I(ia);
        
        I = interp1( theta, I, linspace( -180, 180, 5000 ) )';
        I(find(isnan(I))) = 0;
        
%         I = Gaussian_blurr(I, 100, 10);
        sum_I = sum_I + I;
        
        

        
    end
    
end


end

