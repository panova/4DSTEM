function [ Edge_map ] = Correlation_map_edges( Corr_edges )
%CORRELATION_MAP_EDGES Summary of this function goes here
%   Detailed explanation goes here

Corr_edges(1:5:end, :) = [];
Corr_edges(:, 1:5:end) = [];

Corr_edges(1:4, :) = -1;
Corr_edges(:, 1:4) = -1;
Corr_edges(end-4:end, :) = -1;
Corr_edges(:, end-4:end) = -1;

figure(567); 
clf(); 
imagesc(Corr_edges),
colormap(violetFire(256));
axis equal off; 
caxis([0 1]);

Edge_map = Corr_edges;

find(Corr_edges == 0);
[BlanksRow, BlanksCol] = find(Corr_edges == 0);

for i = 1:numel(BlanksRow)

   left = Corr_edges(BlanksRow(i), BlanksCol(i)-1);
   right = Corr_edges(BlanksRow(i), BlanksCol(i)+1);
   up = Corr_edges(BlanksRow(i)-1, BlanksCol(i));
   down = Corr_edges(BlanksRow(i)+1, BlanksCol(i));
   upleft = Corr_edges(BlanksRow(i)-1, BlanksCol(i)-1);
   upright = Corr_edges(BlanksRow(i)-1, BlanksCol(i)+1);
   downleft = Corr_edges(BlanksRow(i)+1, BlanksCol(i)-1);
   downright = Corr_edges(BlanksRow(i)+1, BlanksCol(i)+1);
   
   neighbors = [left, right, up, down, upleft, upright, downleft, downright];
   neighbors = neighbors(neighbors ~=0);
   
   Edge_map(BlanksRow(i), BlanksCol(i)) = mean(neighbors);


end


figure(567); 
clf(); 
imagesc(1-Edge_map),
colormap(violetFire(256));
% colormap(jet)
axis equal off; 
caxis([0 1]);
end

