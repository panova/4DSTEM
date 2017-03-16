function [  ] = Radon_transform( I )
%RADON_TRANSFORM Summary of this function goes here
%   Detailed explanation goes here

dimension = size(I);
g = fspecial('gaussian', [dimension(1), dimension(2)], dimension(1)/5);
g_in = fspecial('gaussian', [dimension(1), dimension(2)], dimension(1)/5);

g = g/max(g(:));
g_in = g_in/max(g_in(:));

g = g_in.*(1-g);
g = g/max(g(:));
g = g.^5;
I = (I.*g);
theta = 1:360;

[R, xp] = radon(I, theta);
R = R/max(R(:));

figure(77);
clf();
imagesc(theta, xp, R); xlabel('\theta (degrees)');
colormap(jet)

figure(49);
clf();
imagesc(I);
axis equal off;
colormap(jet)


figure(453);
clf();
imagesc(g);
axis equal off;
colormap(jet)



figure(488);
clf();
imagesc(imrotate(I, 180-239));
axis equal off;
colormap(jet)




end

