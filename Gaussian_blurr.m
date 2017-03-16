function [ blured_image ] = Gaussian_blurr( I, size, std )
%Small function that performs a blurr of a 2D image by convolving it with
%a square 2D gaussian of size "size" and standard deviation "std".

gauss = fspecial('gaussian',  [size size],  std);
blured_image = conv2(I, gauss, 'same');


end

