thetas4 = FFATCDIO4(:, :, 1);
thetas22 = FFATC22(:, :, 1);

thetas4 = reshape(thetas4, [1, 128*128*4]);
thetas22 = reshape(thetas22, [1, 128*128*4]);

thetas4 = thetas4(thetas4 ~= -1);
thetas22 = thetas22(thetas22 ~= -1);

h4 = histogram(thetas4, 180);
h22 = histogram(thetas22, 180);

counts4 = h4.Values;
counts22 = h22.Values;

gauss4 = conv(counts4, fspecial('gaussian', [1 5], 1), 'same');
gauss22 = conv(counts22, fspecial('gaussian', [1 5], 1), 'same');

figure(4); clf(); hold on; hist(thetas4, 180); plot(1:180, gauss4, 'r', 'LineWidth', 3); plot(1:179, (diff(gauss4)+200), 'g', 'LineWidth', 3); colormap(jet);
figure(22); clf(); hold on; hist(thetas22, 180); plot(1:180, gauss22,  'r', 'LineWidth', 3); plot(1:179, (diff(gauss22)+200), 'g', 'LineWidth', 3); colormap(jet);


