function [ shift_STACK, distances, del_rows, del_cols ] = Ellipse_realignment_of_stack( STACK )
%ELLIPSE_REALIGNMENT_OF_STACK Summary of this function goes here

dimension = size(STACK);
shift_STACK = zeros([dimension(1), dimension(2), dimension(3)]);

%The center of reference that will become the center of all of the DPs is
%just the geometric middle of the image field of view 
def_center = [round(dimension(1)/2), round(dimension(2)/2)];

%This array stores the distance of the fitted ellipse center to that of the
%reference center: 
distances = zeros([dimension(3), 1]);
del_rows = zeros([dimension(3), 1]);
del_cols = zeros([dimension(3), 1]);


del_row = 0;
del_col = 0;

for ID = 1:dimension(3);% 60*64+50:60*64+50%dimension(3);
   I = STACK(:, :, ID);
   ID

   %The ellipse fitting algorithm, with the row and column coordinates of
   %the fitted ellipse center:
   [ center_row, center_col ] = Ellipse_fitting( I );
   if center_row == 0
%        shift_STACK(:, :, ID) = circshift(I, [del_row, del_col]);
       distances(ID) = sqrt(  (del_row)^2 + (del_col)^2  );
   else 
   %Calculation of the distance of the fitted ellipse center to the
   %reference center:
   del_row = round(def_center(1) - center_row);
   del_col = round(def_center(2) - center_col);
   del_rows(ID) = del_row;
   del_cols(ID) = del_col;

   distances(ID) = sqrt(  (del_row)^2 + (del_col)^2  );
   %shifting the DP: 
%    shift_STACK(:, :, ID) = circshift(I, [del_row, del_col]);
   
%    I = circshift(I, [del_row, del_col]);
%    I(def_center(1), def_center(2)) = 500;
%    figure(12);
%    imagesc(circshift(I, [del_row, del_col]));
%    axis equal off;
   end
end


distances = reshape(distances, [sqrt(dimension(3)), sqrt(dimension(3))]);
del_rows = reshape(del_rows, [sqrt(dimension(3)), sqrt(dimension(3))]);
del_cols = reshape(del_cols, [sqrt(dimension(3)), sqrt(dimension(3))]);
% 
% %We fit the row-shifts and column-shifts to a parabola:
[col, row] = meshgrid(linspace(1, sqrt(dimension(3)), sqrt(dimension(3))), linspace(1, sqrt(dimension(3)), sqrt(dimension(3))));
[col2, row2, del_rows] = prepareCurveData(col, row, medfilt2(del_rows, [5, 5]));
[col2, row2, del_cols] = prepareCurveData(col, row, medfilt2(del_cols, [5, 5]));

fitrowshift = fit([col2,row2],del_rows,'poly22');
fitcolshift = fit([col2,row2],del_cols,'poly22');

fitrowshift = round(fitrowshift([col2, row2]));
fitcolshift = round(fitcolshift([col2, row2]));

for ID = 1:dimension(3);
    ID
    shift_STACK(:, :, ID) = circshift(STACK(:, :, ID), [fitrowshift(ID), fitcolshift(ID)]);
end

fitrowshift = reshape(fitrowshift, [sqrt(dimension(3)), sqrt(dimension(3))]);
fitcolshift = reshape(fitcolshift, [sqrt(dimension(3)), sqrt(dimension(3))]);

del_rows = reshape(del_rows, [sqrt(dimension(3)), sqrt(dimension(3))]);
del_cols = reshape(del_cols, [sqrt(dimension(3)), sqrt(dimension(3))]);


figure(1); clf(); imagesc(medfilt2(del_rows, [5 5])); axis equal off; colormap(jet); caxis([-4 7]);
figure(2); clf(); imagesc(fitrowshift);               axis equal off; colormap(jet); caxis([-4 7]);

figure(11); clf(); imagesc(medfilt2(del_cols, [5 5])); axis equal off; colormap(jet); caxis([-4 7]);
figure(22); clf(); imagesc(fitcolshift);               axis equal off; colormap(jet); caxis([-4 7]);


figure(12);
clf();
imagesc(mean(shift_STACK, 3));
axis equal off;


end

