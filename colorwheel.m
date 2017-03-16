function [  ] = colorwheel( optional_colormap )
%CREATES A COLOWHEEL FOR ANGLES 

colormap default;
jet_wrap = 0;
N = 256;                                % Number Of Colours In Colour Wheel
N = 128;
th = linspace(0, 2*pi, N);              % Angles
xy = @(a,r) [r.*cos(a); r.*sin(a)];     % Function: Calculate (x,y) Vectors
r_ext = 1.0;                            % Outer Radius
r_int = 0.7;                            % Inner Radius
C(1,:,:) = xy(th,r_ext);                % Outer Circle
C(2,:,:) = xy(th,r_int);                % Inner Circle
size(jet);
jet_wrap = vertcat(jet,flipud(jet));
jet_wrap = vertcat(jet_wrap,flipud(jet_wrap));
if nargin >0
    jet_wrap = vertcat(optional_colormap, (optional_colormap));
%     jet_wrap = (vertcat(jet_wrap, flipud(jet_wrap)));
end
size(jet_wrap)
% c = colormap(jet_wrap(N));                   % Set ‘colormap’
figure(2)
C1 = squeeze(C(1,:,:));                 % Reduce ‘C1’ Dimensions
plot(C1(1,:), C1(2,:))
hold on
C2 = squeeze(C(2,:,:));                 % Reduce ‘C2’ Dimensions
plot(C2(1,:), C2(2,:))
for k1 = 1:size(C2,2)
%     plot([C1(1,k1), C2(1,k1)], [C1(2,k1), C2(2,k1)], 'Color',jet_wrap(k1,:), 'LineWidth',0.25)
    k2 = max([mod(k1+1,size(C2,2)), 1]);
    patch([C1(1,k1), C1(1,k2), C2(1,k2), C2(1,k1)], [C1(2,k1), C1(2,k2), C2(2,k2), C2(2,k1)], jet_wrap(k1,:), 'EdgeColor','none')    
end
hold off
axis square off


end

