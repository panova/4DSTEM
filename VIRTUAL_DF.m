function [ DF ] = VIRTUAL_DF( STACK, template, dim1, dim2 )
% This function takes a 4DSTEM stack and performs a simple summation 
% virtual dark field given a template. 
% The convention is such that the DP is hidden where the template is equal 
% to zero and is summed over the locations where the template is = 1 (or
% not equal to zero).
% dim1 and dim2 are the dimensions in real space. 

dimension = size(STACK);

template = repmat(template, [1, 1, dimension(3)]);

% STACK(STACK<0) = 0;

STACK = STACK.*template;

%This is for the format of a three-dimensional stack, where it has not 
% been reshaped into a 4-D stack

DF = sum(sum(STACK, 1), 2);

DF = double(DF);
DF = reshape(DF, [dim1, dim2]);
% DF = DF + abs(min(DF(:)));
DF = DF/abs(max(DF(:)));

figure(4536);
clf();
imagesc(DF)
axis equal off;


end

