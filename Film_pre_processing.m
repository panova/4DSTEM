function [ binary_locations, scan_means, first_frame ] = Film_pre_processing( film, center )

dimension = size(film);

film_cut = film(center(1)-70:center(1)+70, center(2)-70:center(2)+70, :);
film_cut = film;

% 
% figure(3);
% clf();
% imagesc(film_cut(:, :, 10));
% axis equal off

%Taking the sum of all the pixels for each frame:
s = squeeze(sum(sum(film, 1), 2));
s = squeeze(sum(sum(film_cut, 1), 2));
s = s - mean(s(:));

df = 1/dimension(3);
f = -1/2:df:1/2-df;

ff =  abs(fftshift(fft(fftshift(s))));
ff(f<0.01) = 0;
p = find(ff == max(ff(:)));
periodicity = max(1./f(p))

% 
% figure(4);
% clf();
% plot(f, ff);

%Convolving with a gaussian to take away the noise, such that transition
%peaks can be made aparent 
sconv = conv(s, fspecial('gaussian', [round(periodicity/2), 1], 5), 'same');
sconvRIGHT = circshift(sconv, 1);
sconvLEFT = circshift(sconv, -1);

peaks = (((sconv - sconvLEFT)<0)*1).*(((sconv - sconvRIGHT)<0)*1);

%Locations of all the peaks found on the convolved function:
binary_locations = find(peaks);

%Number of peaks found in total. This is the number of transitions, aka
%number of scans. 
peak_number = size(binary_locations);


for scan = 1:peak_number;
    %The location of this scan's peak in the dataset:
    location = binary_locations(scan);
    try
        %Cutting out a window around the peak:
        window = s((location-round(periodicity/2)):(location+round(periodicity/2)));
    catch
        window = s(location);
    end
    
    %location of the true peak within the window. We only want the leftmost
    %peak (to palliate to plateaus), hence the "1".
    true_peak = find(   (window == min(window(:))), 1   );
    delta_location = true_peak - (round(periodicity/2)+1);
    binary_locations(scan) = binary_locations(scan) + delta_location;
end

binary_locations = binary_locations(binary_locations > 0);
binary_locations = unique(binary_locations);
num_scans = size(binary_locations);



%====================================================================
%Adding all of the frames from the same scan together:
adjusted_binary_locations = zeros([num_scans]) + binary_locations;

for scan = 1:num_scans;
    t = circshift(binary_locations, 1);
    diff = abs(t(scan) - binary_locations(scan));
    if diff <= periodicity - 0.5*periodicity
        loc = s(binary_locations(scan));
        locR = s(t(scan));
        if loc > locR
            adjusted_binary_locations(scan) = -15;
        end
        if loc <= locR
            adjusted_binary_locations(scan+1) = -15;
        end
    end
end

adjusted_binary_locations(adjusted_binary_locations == -15) = [];

num_scans = size(adjusted_binary_locations);

readjusted_binary_locations = squeeze(adjusted_binary_locations');
for scan = 2:num_scans-1;
    t = circshift(adjusted_binary_locations, -1);
    diff = abs(t(scan) - adjusted_binary_locations(scan));
    if diff > periodicity + 0.5*periodicity
        left = adjusted_binary_locations(scan)+round(periodicity-periodicity*0.3);
        right = adjusted_binary_locations(scan)+round(periodicity+periodicity*0.3);
        new_loc = find(  (s(left:right) == min(s(left:right))), 1  );
        new_peak = squeeze(new_loc + adjusted_binary_locations(scan)+ round(periodicity-periodicity*0.5));
        readjusted_binary_locations = [readjusted_binary_locations, new_peak];
    end
end

readjusted_binary_locations = sort(readjusted_binary_locations);

bl = abs(readjusted_binary_locations' - circshift(readjusted_binary_locations', 1));

num_scans = size(readjusted_binary_locations');
scan_means = zeros([dimension(1), dimension(2), (num_scans(1))]);

for scan = 1:num_scans-1;
    %This is if we want to average over the WHOLE scan (leads to a lot of
    %spot superimpositions)
    scan_mean = mean(film(:, :, readjusted_binary_locations(scan):readjusted_binary_locations(scan+1)), 3);
    %This is only averaging the first 3 frames of the scan, leading to
    %lower signal but better spot definition/tracking:
%     scan_mean = mean(film(:, :, readjusted_binary_locations(scan):readjusted_binary_locations(scan)+3), 3);

    scan_means(:, :, scan) = scan_mean;
    
end

first_frame = mean(film(:, :, readjusted_binary_locations(1):readjusted_binary_locations(1)), 3);

figure(4);
clf();
imagesc(first_frame);
axis equal off;


scopy = s;
scopy(readjusted_binary_locations) = 1;
scopy(scopy~=1) = 0;
scopy = scopy.*s;

figure(11);
clf();
hold on;
plot(s, 'b-');
% plot(scopy, 'r.');
plot(readjusted_binary_locations, scopy(readjusted_binary_locations), 'r.');
hold off;
xlabel('Frame number');
ylabel('Average intensity over frame');
xlim([0; dimension(3)]);
ylim([min(s(:)), max(s(:))]);

% figure(5);
% clf();
% imagesc(scan_means(:, :, 7));
% axis equal off;


% for i = 1:(num_scans(1)-1); 
%     figure(10);
%     clf();
%     imagesc(scan_means(:, :, i)); 
%     colormap(violetFire(256));
%     title(i);
%     axis equal off; 
%     pause(0);
% end

end

