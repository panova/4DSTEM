function [ Corr_map, Corr_map_edge  ] = DP_cross_correlation( STACK, realDimX, realDimY, center, mycmap )
%GRAIN_SEPARATION_FLOODING Summary of this function goes here
%   Detailed explanation goes here

stack_dimension = size(STACK);
% padd = 3;

lattice_rad = 76;
spot_rad = 12;

% mean_stack = mean(STACK, 3);

STACK = reshape(STACK, [stack_dimension(1), stack_dimension(2), realDimX, realDimY]);
% STACK = padarray(STACK, [0, 0, padd, padd], 'both');

stack_dimension = size(STACK);

Corr_map = zeros([realDimX.*3, realDimY.*3]);
Corr_map_edge = zeros([realDimX.*5, realDimY.*5]);

[Ya, Xa] = meshgrid(1:stack_dimension(2), 1:stack_dimension(1));
Ra = sqrt((center(1) - Xa).^2 + (center(2)-Ya).^2);
mask = (lattice_rad) - Ra;
mask2 = (lattice_rad + 3*spot_rad) - Ra;
mask(mask>=0) = 1;
mask(mask<0) = 0;
mask2(mask2>=0) = 1;
mask2(mask2<0) = 0;
aperture = abs((mask - mask2));

tic
for R = 2:(stack_dimension(3)-1)
    
    R
    for C = 2:(stack_dimension(4)-1)     
        
        DP = STACK(:, :, R, C);
        newR = (R.*3) - 1;
        newC = (C.*3) - 1;
        
        newR_edge = (R.*5) - 2;
        newC_edge = (C.*5) - 2;        
        
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
        
        I = interp1( theta, I, linspace( -180, 180, numel(I) ) )';
        I(find(isnan(I))) = 0;

        I = Gaussian_blurr(I, 100, 20);

                
%         figure(342)
%         clf();
%         plot(theta, I, '-');
%         title('I at R,C');
%         drawnow;
%         
        
        %Direct neighbors:
        %LEFT
            pleft = STACK(:, :, R, C-1);
%             pleft = mean_stack;
            %masking anything that is NOT the annulus whithin which
            %diffraction is contained:
            pleft = pleft.*aperture;
            
            %Locating all the points within the image that have not been
            %erased by the aperture:
            [rowsleft, colsleft, Ileft] = find(pleft);
            
            %Shifting them to the center of the image so that the angle can
            %be accurately computed for each point:
            rowsleft = rowsleft - center(1);
            colsleft = colsleft - center(2);
            
            %Computing the angle theta at which each point is located:
            thetaleft = rad2deg(atan2(rowsleft, colsleft));
            
            %The points are not in order of ascending angle. We make it so:
            ptsleft = sortrows(cat(2, thetaleft, Ileft));  
            
            %Now the points are in order, but there are duplicates. We
            %suppress them:
            Ileft = ptsleft(:, 2);
            [thetaleft, ialeft, ~] = unique(ptsleft(:, 1), 'stable');
            Ileft = Ileft(ialeft);
                        
            %Because the number of points included within our aperture is
            %arbitrary, we want the two signal curves we are correlating to
            %be the same size. However, we can't just truncate it - that
            %would be random and the two curves would not correspond to one
            %another. So we want to re-sample the second curve over the
            %range of angles (-180 to 180) with the same number points as
            %our base curve, I. Now that I think of it, for consistency's
            %sake, one should ALSO resample the base curve, to have equally
            %separated points (because right now they are plotted as a
            %function of the calculated angle at which the point is, which
            %absolutely does not need to be equally spaced). So I am doing
            %that above too, but this is the explanation of why. 
            I_interpleft = interp1( thetaleft, Ileft, linspace( -180, 180, numel(I) ) )';
            %So apparently interp1 gets all f*cky when it tries to
            %interpolate data outside of the domain of the given function it
            %tries to fit, so, well, it sticks a bunch of "NaNs" in
            %there at the tails. Because you can't divide by NaN without
            %breaking the universe, you have to suppress them, so here we
            %go:
            I_interpleft(find(isnan(I_interpleft))) = 0;
            
            %Finally we blurr the hell out if it because who needs noisy
            %data:
            I_interpleft = Gaussian_blurr(I_interpleft, 100, 20);

            
            [Correlation_coeff0, ~] = corrcoef(I, I_interpleft);    
            Corr_map(newR, newC-1) = Correlation_coeff0(1, 2);
            Corr_map_edge(newR_edge, newC_edge-1) = Correlation_coeff0(1, 2);
            
%             
%             figure(4567);
%             clf();
%             plot(linspace(-180, 180, numel(I)), I_interpleft);
%             title('I interpret LEFT'); 
%             
%             
%             figure(57);
%             clf();
%             plot(xcorr(I, I_interpleft));
%             title('correlation');
            
        
        %RIGHT
            pright = STACK(:, :, R, C+1);
            pright = pright.*aperture;
            [rowsright, colsright, Iright] = find(pright);
            rowsright = rowsright - center(1);
            colsright = colsright - center(2);
            thetaright = rad2deg(atan2(rowsright, colsright));
            ptsright = sortrows(cat(2, thetaright, Iright));
            Iright = ptsright(:, 2);
            [thetaright, iaright, ~] = unique(ptsright(:, 1), 'stable');
            Iright = Iright(iaright);
            I_interpright = interp1( thetaright, Iright, linspace( -180, 180, numel(I) ) )';
            I_interpright(find(isnan(I_interpright))) = 0;
            I_interpright = Gaussian_blurr(I_interpright, 100, 20);
            [Correlation_coeff1, ~] = corrcoef(I, I_interpright);
            Corr_map(newR, newC+1) = Correlation_coeff1(1, 2);
            Corr_map_edge(newR_edge, newC_edge+1) = Correlation_coeff1(1, 2);

        
        %UP
            %So at this point I got really super bored of replacing all the
            %"lefts" with "rights" and "ups" so since all of these are
            %getting overwritten anyway, I'm just keeping all the variables as "right" 
            pup = STACK(:, :, R-1, C);
            pright = pup.*aperture;
            [rowsright, colsright, Iright] = find(pright);
            rowsright = rowsright - center(1);
            colsright = colsright - center(2);
            thetaright = rad2deg(atan2(rowsright, colsright));
            ptsright = sortrows(cat(2, thetaright, Iright));
            Iright = ptsright(:, 2);
            [thetaright, iaright, ~] = unique(ptsright(:, 1), 'stable');
            Iright = Iright(iaright);
            I_interpright = interp1( thetaright, Iright, linspace( -180, 180, numel(I) ) )';
            I_interpright(find(isnan(I_interpright))) = 0;
            I_interpright = Gaussian_blurr(I_interpright, 100, 20);
            [Correlation_coeff2, ~] = corrcoef(I, I_interpright);
            Corr_map(newR-1, newC) = Correlation_coeff2(1, 2);
            Corr_map_edge(newR_edge-1, newC_edge) = Correlation_coeff2(1, 2);

        
        %DOWN
            pdown = STACK(:, :, R+1, C);
            pright = pdown.*aperture;
            [rowsright, colsright, Iright] = find(pright);
            rowsright = rowsright - center(1);
            colsright = colsright - center(2);
            thetaright = rad2deg(atan2(rowsright, colsright));
            ptsright = sortrows(cat(2, thetaright, Iright));
            Iright = ptsright(:, 2);
            [thetaright, iaright, ~] = unique(ptsright(:, 1), 'stable');
            Iright = Iright(iaright);
            I_interpright = interp1( thetaright, Iright, linspace( -180, 180, numel(I) ) )';
            I_interpright(find(isnan(I_interpright))) = 0;
            I_interpright = Gaussian_blurr(I_interpright, 100, 20);
            [Correlation_coeff3, ~] = corrcoef(I, I_interpright);
            Corr_map(newR+1, newC) = Correlation_coeff3(1, 2); 
            Corr_map_edge(newR_edge+1, newC_edge) = Correlation_coeff3(1, 2);
            
        %Diagonal neighbors:
        
        %LEFT UP
            pleftup = STACK(:, :, R-1, C-1);
            pright = pleftup.*aperture;
            [rowsright, colsright, Iright] = find(pright);
            rowsright = rowsright - center(1);
            colsright = colsright - center(2);
            thetaright = rad2deg(atan2(rowsright, colsright));
            ptsright = sortrows(cat(2, thetaright, Iright));
            Iright = ptsright(:, 2);
            [thetaright, iaright, ~] = unique(ptsright(:, 1), 'stable');
            Iright = Iright(iaright);
            I_interpright = interp1( thetaright, Iright, linspace( -180, 180, numel(I) ) )';
            I_interpright(find(isnan(I_interpright))) = 0;
            I_interpright = Gaussian_blurr(I_interpright, 100, 20);
            [Correlation_coeff4, ~] = corrcoef(I, I_interpright);
            Corr_map(newR-1, newC-1) = Correlation_coeff4(1, 2);
            Corr_map_edge(newR_edge-1, newC_edge-1) = Correlation_coeff4(1, 2);
            
        %RIGHT UP
            prightup = STACK(:, :, R-1, C+1);
            pright = prightup.*aperture;
            [rowsright, colsright, Iright] = find(pright);
            rowsright = rowsright - center(1);
            colsright = colsright - center(2);
            thetaright = rad2deg(atan2(rowsright, colsright));
            ptsright = sortrows(cat(2, thetaright, Iright));
            Iright = ptsright(:, 2);
            [thetaright, iaright, ~] = unique(ptsright(:, 1), 'stable');
            Iright = Iright(iaright);
            I_interpright = interp1( thetaright, Iright, linspace( -180, 180, numel(I) ) )';
            I_interpright(find(isnan(I_interpright))) = 0;
            I_interpright = Gaussian_blurr(I_interpright, 100, 20);
            [Correlation_coeff5, ~] = corrcoef(I, I_interpright);
            Corr_map(newR-1, newC+1) = Correlation_coeff5(1, 2);
            Corr_map_edge(newR_edge-1, newC_edge+1) = Correlation_coeff5(1, 2);
        
        %LEFT DOWN
            pleftdown = STACK(:, :, R+1, C-1);
            pright = pleftdown.*aperture;
            [rowsright, colsright, Iright] = find(pright);
            rowsright = rowsright - center(1);
            colsright = colsright - center(2);
            thetaright = rad2deg(atan2(rowsright, colsright));
            ptsright = sortrows(cat(2, thetaright, Iright));
            Iright = ptsright(:, 2);
            [thetaright, iaright, ~] = unique(ptsright(:, 1), 'stable');
            Iright = Iright(iaright);
            I_interpright = interp1( thetaright, Iright, linspace( -180, 180, numel(I) ) )';
            I_interpright(find(isnan(I_interpright))) = 0;
            I_interpright = Gaussian_blurr(I_interpright, 100, 20);
            [Correlation_coeff6, ~] = corrcoef(I, I_interpright);
            Corr_map(newR+1, newC-1) = Correlation_coeff6(1, 2);
            Corr_map_edge(newR_edge+1, newC_edge-1) = Correlation_coeff6(1, 2);
        

        %RIGHT DOWN
            prightdown = STACK(:, :, R+1, C+1);
            pright = prightdown.*aperture;
            [rowsright, colsright, Iright] = find(pright);
            rowsright = rowsright - center(1);
            colsright = colsright - center(2);
            thetaright = rad2deg(atan2(rowsright, colsright));
            ptsright = sortrows(cat(2, thetaright, Iright));
            Iright = ptsright(:, 2);
            [thetaright, iaright, ~] = unique(ptsright(:, 1), 'stable');
            Iright = Iright(iaright);
            I_interpright = interp1( thetaright, Iright, linspace( -180, 180, numel(I) ) )';
            I_interpright(find(isnan(I_interpright))) = 0;
            I_interpright = Gaussian_blurr(I_interpright, 100, 20);
            [Correlation_coeff7, ~] = corrcoef(I, I_interpright);
            Corr_map(newR+1, newC+1) = Correlation_coeff7(1, 2);
            Corr_map_edge(newR_edge+1, newC_edge+1) = Correlation_coeff7(1, 2);
            
            Corr_map(newR, newC) = mean([Correlation_coeff0(1, 2), Correlation_coeff1(1, 2), Correlation_coeff2(1, 2), Correlation_coeff3(1, 2), Correlation_coeff4(1, 2), Correlation_coeff5(1, 2), Correlation_coeff6(1, 2), Correlation_coeff7(1, 2)]);
            Corr_map_edge(newR_edge, newC_edge) = mean([Correlation_coeff0(1, 2), Correlation_coeff1(1, 2), Correlation_coeff2(1, 2), Correlation_coeff3(1, 2), Correlation_coeff4(1, 2), Correlation_coeff5(1, 2), Correlation_coeff6(1, 2), Correlation_coeff7(1, 2)]);


    end
    
end

toc

figure(101);
clf();
imagesc(Corr_map);
%caxis([min(Grain_map(Grain_map>0)), max(Grain_map(Grain_map>0))]);
% caxis([0 180])
title('Correlation Map');
colormap(violetFire(256));
axis equal off
drawnow;


figure(102);
clf();
imagesc(Corr_map_edge);
%caxis([min(Grain_map(Grain_map>0)), max(Grain_map(Grain_map>0))]);
% caxis([0 180])
title('Correlation Map with Edges');
colormap(violetFire(256));
axis equal off
drawnow;




end



