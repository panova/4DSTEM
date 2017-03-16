function [  ] = PolyFit2D( I )
%POLYFIT2D Summary of this function goes here
%   Detailed explanation goes here

dimension = size(I);

[x, y] = meshgrid(1:dimension(1), 1:dimension(2));

x = reshape(x, [1, dimension(1)*dimension(2)]);
y = reshape(y, [1, dimension(1)*dimension(2)]);
I = reshape(I, [1, dimension(1)*dimension(2)]);

[x2, y2, I2] = prepareSurfaceData(x, y, I);

fitsurface = fit([x2,y2],I2, 'poly22','Normalize','on');

figure(2)
clf();
plot(fitsurface)





end

