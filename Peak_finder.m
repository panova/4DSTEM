function [ peak_list ] = Peak_finder( I, number_of_peaks_to_return )
%Peak search (searching for pixels brighter than their surrounding
%pixels).

%INPUT:   -- I (we want this function to find peaks on a 2D image regardless of
%what that image is); this function only finds peaks - it starts with the
%assumption that all the necessary blurring/cleaning/fftshift operations have
%already been performed. 

%RETURNS: -- List of peak positions along with the intensity values of
%those peaks

%Note: the original PC can be easily accessed via
%[PC, template] = PhaseCorellation(I);

%Eliminating small noisy peaks:
% I(I < max(I(:))*0.01) = 0;
dimension = size(I);
center = [round(dimension(1)/2) round(dimension(2)/2)];

%Adding a small amount of noise to the image so that plateaus (peaks that
%span several pixels) can be found (plateaus will not be found via 
%the technique developped below).
diff_noise = randn(dimension(1), dimension(2))*max(I(:))/10000;
I = I + diff_noise;

%Shifted versions of PC for comparison of neighbor values:
left = circshift(I, [0, -1]);
right = circshift(I, [0, 1]);
top = circshift(I, [-1, 0]);
bottom = circshift(I, [1, 0]);

%Comparing elements of PC to their four nearest neighbors:
peak = ((I > left) & (I > right) & (I > top) & (I > bottom))*1;
%Masking the original PC array with the boolean "peak" array obtained
%above:
I_peak = I;
I_peak(peak==0) = 0;

peak_list = zeros(number_of_peaks_to_return, 3);

%If only one peak is wanted per image, the strongest peak is stored along 
%with its value on the image:
if number_of_peaks_to_return == 1
    [peak_row, peak_col] = find(I_peak==max(I_peak(:)));
    if size(peak_row, 1) > 1
        peak_row = center(1);
        %disp('oops');
    end
    if size(peak_col, 1) > 1
        peak_col = center(2);
    end
    peak_list = [peak_row, peak_col, max(I_peak(:))];
end


%If several peaks are needed (more than one disk per image):
if number_of_peaks_to_return ~=1
    for step = 1:number_of_peaks_to_return
        [peak_row, peak_col] = find(I_peak==max(I_peak(:)));
        peak_list(step, :) = [peak_row, peak_col, max(I_peak(:))];
        I_peak(peak_row, peak_col) = 0;
    end
end


end

