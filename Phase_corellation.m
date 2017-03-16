function [PC, template] = Phase_corellation( I, rad, temp )
%This function computes the phase correlation of an image with a template and uses it to
%find disks on the image.

dimension = size(I);
%radius of the template disk, measured approximatively by eye
temp_radius = rad;
%position of the center of the template disk
center = [round(dimension(1)/2) round(dimension(2)/2)];

%Tracing the circle onto the template:
[Ya, Xa] = meshgrid(1:dimension(2), 1:dimension(1));
Ra = sqrt((center(1) - Xa).^2 + (center(2)-Ya).^2);
template = temp_radius - Ra + 0.5;
template(template<0) = 0;
template(template>1) = 1;

if nargin > 2;
    template = temp;
end

template = ifftshift(template);

%The normalization of the phase correlation:
norml = abs(fft2(I).*conj(fft2(template)));
%Computing the actual phase correlation:
PC = abs(ifft2(fft2(I).*conj(fft2(template))./norml));

% figure(37)
% clf
% imagesc(PC);
% axis equal off
% colormap(violetFire(256));


end


