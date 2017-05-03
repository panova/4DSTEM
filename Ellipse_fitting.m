function [ center_row, center_col ] = Ellipse_fitting( I )
%This function takes in a DP as input and fits an ellipse to its 
%halo in order to determine the center of the beam. We assume that the
%images are integers, varying from values 1 -> 255. 

dimension = size(I);

Im = I;
I = single(I);

I = I/255;

% I = I/max(I(:));


% Thresholding the image to get only a partial ring from the halo
% for E1:
% I = (abs(I/max(I(:)) - 0.3)<0.025)*1;
% I = (abs(I/max(I(:)) - 0.2)<0.025)*1;

%for D2:
% g = fspecial('gaussian', [dimension(1), dimension(2)], dimension(1)/5);
% g = g/max(g(:));
% I = I.*g;
% I = (abs(I/max(I(:)) - 0.2)<0.01)*1;

%FOR FILMS PE 2016:
%film10
% g = fspecial('gaussian', [dimension(1), dimension(2)], dimension(1)/10);
% g = circshift(g, [65 0]);
% g = g/max(g(:));
% I = I.*g;
% I = I./max(I(:));

%For TACAKS TITAN X datasets, the sigma for the gaussian is dimension(1)/15
%For Xi datasets, the sigma for the gaussian is dimension(1)/10
g = fspecial('gaussian', [dimension(1), dimension(2)], dimension(1)/10);
g = g/max(g(:));
I = I.*g;
% I = I./max(I(:));

% figure(33);
% imagesc(I);
% colormap(gray);
% axis equal off;
% drawnow;

%To figure out the numbers:
%Plot I, the image you want to threshold. Find the min and the max values
%that you want to include within your threshold. The difference in those
%values is what you want to subtract from I. 
%Then take the absolute value. 
%Then set the threshold limit 
%For TACAKS TITAN:
% I = ( abs(I-(0.15-0.04)) <0.02 )*1;

%For Xi15: 
I = ( abs(I-(0.3-0.01)) <0.02 )*1;

% I = (     abs(I - 0.5) < 0.4 )*1;

% I = (     abs(I - 0.7) < 0.2     )*1;
% 
% figure(30);
% imagesc(I);
% colormap(gray);
% axis equal off;
% % caxis([0 0.5])
% drawnow

if sum(I(:)) ~=0;
    
    [row_pts, col_pts] = find(I);
    XY = [row_pts  col_pts];
    ellipse = EllipseDirectFit(XY);
    
    a = ellipse(1);
    b = ellipse(2)/2;
    c = ellipse(3);
    d = ellipse(4)/2;
    e = ellipse(5)/2;
    g = ellipse(6);
    
    center_row = (c*d - b*e) / (b^2 - a*c);
    center_col = (a*e - b*d) / (b^2 - a*c);
    
    Im(round(center_row), round(center_col)) = 500;
    Im(140, 140) = 350;
else
    center_row = 0;
    center_col = 0;
end





%Convert the A to str 
a = num2str(a); 
b = num2str(b*2); 
c = num2str(c); 
d = num2str(d*2); 
e = num2str(e*2); 
g = num2str(g);

eqt= ['(',a, ')*y^2 + (',b,')*x*y + (',c,')*x^2 + (',d,')*y+ (',e,')*x + (',g,')']; 
xmin=0.7*min(XY(:,1)); 
xmax=1.3*max(XY(:,2)); 
% 
% figure(47); 
% clf();
% hold on;
% imagesc(Im);
% axis equal off;
% ex = ezplot(eqt,[xmin,xmax]);
% set(ex, 'color', 'yellow');
% plot(center_col, center_row, 'y+');
% colormap(gray);
% hold off;
% drawnow;

% figure(4); 
% clf();
% hold on;
% imagesc(Im);
% axis equal off;
% scatter(XY(:,2),XY(:,1), 'b.') ;
% ex = ezplot(eqt,[xmin,xmax]);
% set(ex, 'color', 'yellow');
% plot(center_col, center_row, 'y+');
% axis equal off;
% colormap(gray);
% hold off;
% drawnow;


end

