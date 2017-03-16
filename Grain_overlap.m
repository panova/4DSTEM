function [ dimension ] = Grain_overlap( FF1, FF2 )
%GRAIN_OVERLAP Summary of this function goes here
%   Detailed explanation goes here

dimension1 = size(FF1);
dimension2 = size(FF2);

rearrangedFF1 = zeros(size(FF1));
rearrangedFF2 = zeros(size(FF2));

for i = 1:dimension(1)
    for k = 1:dimension(2)
        angle1 = FF1(i, k);
        angle2 = FF2(i, k);
        surr1 = 
    end
end


end

