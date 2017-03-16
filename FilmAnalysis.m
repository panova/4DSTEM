function [ Peaks_ref, Peaks_all, tcs ] = FilmAnalysis( film, start_frame, num_peaks, center, periodicity, PCyn, first_frame )
%FILMANALYSIS Summary of this function goes here
%   Detailed explanation goes here


backup_dim = size(film);
film = film(:, :, start_frame:end-1);
film = 100*film./225;

dimension = size(film);

lattice_rad = 85;

center_shift = [round(dimension(1)/2), round(dimension(2)/2)] - center;

aesthetic_aperture_rad = 14; % This is for display purposes and testing. This is the radius of the expected diffraction spot. 
masking_aperture_rad = aesthetic_aperture_rad*3; %When looking for peaks on the DP, the peaks found will be masked with an aperture of this radius

[column, row] = meshgrid(1:dimension(2), 1:dimension(1));

Peaks_ref = zeros([num_peaks 3]);
Peaks_end_ref = zeros([1, num_peaks]);

%=========================================================================
%SETTING A REFERENCE USING THE FIRST FRAME
ID = 1;
I = film(:, :, ID);

%The first frame may be provided if the average has too many spots or
%washes out the spots too much. The "first frame" argument is usually
%provided by the Film_preprocessing function, but can be set manually. 
if nargin() == 7;
    I = first_frame;
end

I_end = film(:, :, end);
% I = 100*I/255;
Im = I;
% Im = Im/255;
Ir = Im;

Iorig = I;


I = medfilt2(I);
[PC, ~] = Phase_corellation( I, aesthetic_aperture_rad );




%Central aperture masking the bright central beam:
Ra = sqrt((center(1) - row).^2 + (center(2)-column).^2);
central_masked = (lattice_rad/1.3) - Ra + 0.5;
I(central_masked>0) = 0;
PC(central_masked>0) = 0;

Ir(central_masked-20>0) = 0;
%Im is the original image with the center masked:
Im(central_masked>0) = 0;

I(central_masked<-aesthetic_aperture_rad-45) = 0;
PC(central_masked<-aesthetic_aperture_rad-45) = 0;

% Ir(central_masked<-aesthetic_aperture_rad-100) = 0;
Ir(central_masked-20<-aesthetic_aperture_rad) = 0;


k = sum(I(central_masked<-aesthetic_aperture_rad-45));

%Gaussian blurr of the original image:
% I = Gaussian_blurr( I, 10, 3 );

%Finding the peaks on the image:
Ipeaks = zeros(num_peaks, 3);
        
figure(1)
clf()
imagesc(I);
axis equal off

% 
% figure(71);
% clf();
% imagesc((PC));
% axis equal off;


%FINDING PEAKS AND SAVING THEIR LOCATIONS:
for i = 1:num_peaks
       %Actually find the peak
       
       if strcmp(PCyn, 'True') == 1;
           Ipeak  = Peak_finder( PC, 1 );
       else
           Ipeak  = Peak_finder( Gaussian_blurr( I, 10, 3 ), 1 );

       end

       %And then MASK it
       aperture_location = [Ipeak(1, 1) Ipeak(1, 2)];
       Ra = sqrt((aperture_location(1) - row).^2 + (aperture_location(2)-column).^2);
       masked = masking_aperture_rad - Ra + 0.5;
       I(masked>0) = 0;
       PC(masked>0) = 0;
        
       aesthetic = aesthetic_aperture_rad - Ra + 0.5;

       peak_signal = mean(Im(aesthetic>0));%+ sum(Im(aesthetic2>0))
       peak_signal_end = mean(I_end(aesthetic>0));

       Im(aesthetic>0) =  Im(aesthetic>0)+ 50;
       Ir(aesthetic>0) =  0;
% 
%         figure(345);
%         clf();
%         imagesc(Im);
%         axis equal off;   
        
        %re-nolmalizing after the peak has been masked:
        I = I/max(I(:));
        
        Ipeaks(i, 1:2) = Ipeak(1:2);
        Ipeaks(i, 3) = peak_signal;
        
        Peaks_ref(i, 1:2) = Ipeaks(i, 1:2);
        Peaks_ref(i, 3) = Ipeaks(i, 3);
        Peaks_end_ref(i) = peak_signal_end;
end


figure(345);
clf();
imagesc(Im);
axis equal off;

%==========================================================================

figure(12);
clf();
imagesc(Iorig);
hold on;
for i = 1:num_peaks
    text(Peaks_ref(i, 2), Peaks_ref(i, 1), num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
hold off;
axis equal off



% %ROTATING THE FIRST FRAME
% %This is the image that we are rotating, with the peaks masked: 
% Ir = circshift(Ir, [center_shift, center_shift]);
% % pad the image with zeros so we don't lose anything when we rotate.
% [iLength, iWidth] = size(Ir);
% iDiag = sqrt(iLength^2 + iWidth^2);
% LengthPad = ceil(iDiag - iLength) + 2;
% WidthPad = ceil(iDiag - iWidth) + 2;
% padIMG = zeros(iLength+LengthPad, iWidth+WidthPad);
% padIMG(ceil(LengthPad/2):(ceil(LengthPad/2)+iLength-1), ...
%     ceil(WidthPad/2):(ceil(WidthPad/2)+iWidth-1)) = Ir;
%    
% n = size(padIMG,1);
% x = linspace(-1,1,n);
% [X1,Y1] = meshgrid(x,x);
% thetas = 1:360;
% 
% threads = [];
% 
% n = length(thetas);
% PR = zeros(size(padIMG,2), n);
%    
% %For each angle, ROTATE:
% for z = 1:n
%     theta = (90-thetas(z))*pi/180;
%     X = cos(theta)*X1 + -sin(theta)*Y1;
%     Y = sin(theta)*X1 + cos(theta)*Y1;
%     % interpolate
%     tmpimg = interp2(X1,Y1,padIMG,X,Y);
%     tmpimg(isnan(tmpimg)) = 0;
%     
%     thread = tmpimg(:, round(size(tmpimg, 1)/2));
%     
%     %There seems to be TWO characteristic lattice spacings:
%     centval1 = round(size(tmpimg, 1)/2) - lattice_rad-10;
%     centval2 = round(size(tmpimg, 1)/2) - lattice_rad;
%     
%     val1 = mean(thread(centval1-1:centval1+1));
%     val2 = mean(thread(centval2-1:centval2+1));
%     
%     thread = thread(round(size(tmpimg, 1)/2)-lattice_rad-98:round(size(tmpimg, 1)/2));
%     
%     tmpimg(1:round(size(tmpimg, 1)/2), round(size(tmpimg, 1)/2)-3:round(size(tmpimg, 1)/2)+3) = 0;
%     tmpimg(centval1-3:centval1+3, round(size(tmpimg, 1)/2)-3:round(size(tmpimg, 1)/2)+3) = 0.75;
%     tmpimg(centval2-3:centval2+3, round(size(tmpimg, 1)/2)-3:round(size(tmpimg, 1)/2)+3) = 1;
%     
%     %        if val1 == 0 | val2 ==0
%     %
%     %                        figure(45);
%     %                        clf();
%     %                        imagesc(tmpimg);
%     %                        title('REJECTED')
%     %                        axis equal off;
%     %        else
%     %                        figure(47);
%     %                        clf();
%     %                        imagesc(tmpimg);
%     %                        title('accepted')
%     %                        axis equal off;
%     
%     %appending to the threads array the newly created thread:
%     threads = [threads, thread];
%     
% end
%           
% toc
% 
% threads_start = threads;
% size(threads)
% figure(34);
% clf();
% plot(mean(threads(1:140, :), 2), 'r');


%=========================================================================
%Now we plot the statistics of the spots and background:

%Stores the X, Y, signal for each peak; the last peak is the amorphous:
Peaks_all = zeros([num_peaks+1 3 dimension(3)]);

tic
for tempID = 1:dimension(3);
   
   ID = tempID; 
   I = film(:, :, tempID);
%    Im = film(:, :, ID);  
%    I = 100*I/255;

   Im = I;
%    Im = Im/255;
   
   Ir = Im;
   
   %     Central aperture masking the bright central beam:
   Ra = sqrt((center(1) - row).^2 + (center(2)-column).^2);
   central_masked = (lattice_rad/1.3) - Ra + 0.5;
   %Masking the center beam:
   I(central_masked>0) = 0;
   Ir(central_masked-12>0) = 0;

   %Masking the outside of the diffraction ring too, making an annular aperture:
%    I(central_masked<-aesthetic_aperture_rad-45) = 0;
%    Ir(central_masked<-aesthetic_aperture_rad-100) = 0;
   Ir(central_masked-7<-aesthetic_aperture_rad) = 0;

   I(central_masked<-aesthetic_aperture_rad-100) = 0;

%    Ipeaks = zeros(num_peaks, 3);
        
    %Slicing in the Peaks_all array to make parfor happy:
    Peaks_alltemp = Peaks_all(:, :, tempID);
    
    xunit = [];
    yunit = [];
    for i = 1:num_peaks
        Ipeak  = Peaks_ref(i, 1:2);
        aperture_location = [Ipeak(1, 1) Ipeak(1, 2)];
        Ra = sqrt((aperture_location(1) - row).^2 + (aperture_location(2)-column).^2);
        aesthetic = aesthetic_aperture_rad - Ra + 0.5;
%         peak_signal = mean(I(aesthetic>0))/Peaks_ref(i, 3);
        peak_signal = mean(I(aesthetic>0));

        Peaks_alltemp(i, 1:2) = Ipeak;
        Peaks_alltemp(i, 3) = peak_signal;

%         Peaks_all(i, 1:2, ID) = Ipeak;
%         Peaks_all(i, 3, ID) = peak_signal;
        
        I(aesthetic>0) = 0;   
%         Im(aesthetic>0) = 0;
        Im(aesthetic>0) =  Im(aesthetic>0)+ 0;

%         Ir(aesthetic>0) =  0;

        th = 0:pi/50:2*pi;
        xunitf = aesthetic_aperture_rad * cos(th) + Ipeak(1, 2) ;
        yunitf = aesthetic_aperture_rad * sin(th) + Ipeak(1, 1);
        
        
        xunit = cat(1, xunit, xunitf);
        yunit = cat(1, yunit, yunitf);
        
        size(xunit)

    end
        
    %The last entry in Peaks_all is reserved for the amorphous halo value;
    %here we set its "location" to zero (since it's not a peak, there is no
    %x or y to store) and the "amor" (=sum over amorphous) 
    %value is stored in the "intensity" slot
    amor = mean(Ir(Ir>0));
    
    % r= desired radius
    % x = x coordinates of the centroid
    % y = y coordinates of the centroid


    figure(56);
    clf();
    imagesc(Im);
    hold on;
    for i = 1:num_peaks
        plot(xunit(i, :), yunit(i, :), 'r--');
    end
    axis equal off;
    drawnow
    hold off;

    
    Peaks_alltemp(num_peaks+1, 1:2) = [0 0];
    Peaks_alltemp(num_peaks+1, 3) = amor;
    %Here we normalize the peak value (average over spot area) to the
    %average background (average over locations witout spots).Because the
    %background changes as the scans progress, the ratio of spot-to-halo
    %changes too; just normalizing it to the halo on the last frame is not
    %"true". 
    Peaks_alltemp(1:end-1, 3) = Peaks_alltemp(1:end-1, 3)./amor;
    
    figure(390);
    clf();
    hold on
    for i = 1:num_peaks
        plot(linspace(0, (dimension(3).*periodicity/30), dimension(3))', squeeze(Peaks_all(i, 3, :)), 'DisplayName', num2str(i));
    end
    drawnow;
    legend('show');
    hold off;

    
    %Storing the peaks in a way that parfor is happy with:
    Peaks_all(:, :, tempID) = Peaks_alltemp;
    

    
%     if ID == 120
%         tic
%         
%         Ir = circshift(Ir, [center_shift, center_shift]);
%         % pad the image with zeros so we don't lose anything when we rotate.
%         [iLength, iWidth] = size(Ir);
%         iDiag = sqrt(iLength^2 + iWidth^2);
%         LengthPad = ceil(iDiag - iLength) + 2;
%         WidthPad = ceil(iDiag - iWidth) + 2;
%         padIMG = zeros(iLength+LengthPad, iWidth+WidthPad);
%         padIMG(ceil(LengthPad/2):(ceil(LengthPad/2)+iLength-1), ...
%             ceil(WidthPad/2):(ceil(WidthPad/2)+iWidth-1)) = Ir;
%         
%         n = size(padIMG,1);
%         x = linspace(-1,1,n);
%         [X1,Y1] = meshgrid(x,x);
%         thetas = 1:360;
%         threads = [];
%         
%         n = length(thetas);
%         PR = zeros(size(padIMG,2), n);
%         
%         
%         
%         for z = 1:n
%             theta = (90-thetas(z))*pi/180;
%             X = cos(theta)*X1 + -sin(theta)*Y1;
%             Y = sin(theta)*X1 + cos(theta)*Y1;
%             % interpolate
%             tmpimg = interp2(X1,Y1,padIMG,X,Y);
%             tmpimg(isnan(tmpimg)) = 0;
%             
%             thread = tmpimg(:, round(size(tmpimg, 1)/2));
%         
%             %There seems to be TWO characteristic lattice spacings:
%             centval1 = round(size(tmpimg, 1)/2) - lattice_rad-10;
%             centval2 = round(size(tmpimg, 1)/2) - lattice_rad;
%             
%             val1 = mean(thread(centval1-1:centval1+1));
%             val2 = mean(thread(centval2-1:centval2+1));
%             
%             thread = thread(round(size(tmpimg, 1)/2)-lattice_rad-98:round(size(tmpimg, 1)/2));
%             
%             tmpimg(1:round(size(tmpimg, 1)/2), round(size(tmpimg, 1)/2)-3:round(size(tmpimg, 1)/2)+3) = 0;
%             tmpimg(centval1-3:centval1+3, round(size(tmpimg, 1)/2)-3:round(size(tmpimg, 1)/2)+3) = 0.75;
%             tmpimg(centval2-3:centval2+3, round(size(tmpimg, 1)/2)-3:round(size(tmpimg, 1)/2)+3) = 1;
%             
%             if val1 == 0 | val2 ==0
%             
%             else            
%                 threads = [threads, thread];
%             
%             end       
%         end
%         toc
%         figure(37);
%         clf();
%         hold on;
%         plot(mean(threads(2:136, :), 2), 'r');
%         plot(mean(threads_start(2:136, :), 2), 'b');    
% %         fill(mean(threads(2:136, :), 2), mean(threads_start(2:136, :), 2), 'b', 'facealpha', 0.5);
%         hold off;
%         
%     end
    
end

toc

colors = [[0/255,0/255,0/255]; [0/255, 206/255, 209/255]; [25/255, 25/255, 112/255]; [220/255, 20/255, 60/255]; [255/255, 127/255, 80/255]; [128/255, 0/255, 128/255]; [255/255, 20/255, 147/255]];

x = (1:dimension(3))';
secds = linspace(0, (dimension(3).*periodicity/30), dimension(3))' + start_frame.*periodicity/30;
secdsback = linspace(0, (dimension(3).*periodicity/30), dimension(3))';

figure(37);
clf();
col = squeeze(datasample(colors, num_peaks, 'Replace', false)');

tcs = zeros([3, num_peaks]);


for i = 1:num_peaks;
    hold on;
    try
        tofit = squeeze(Peaks_all(i, 3, :))-Peaks_all(i, 3, end);
        
        %For plotting purposes, we keep the whole array:
        tofitfull = tofit;
        tofitfull = tofitfull./max(tofitfull(:));

%         tofit = conv(tofit, fspecial('gaussian', [5,1], 3), 'same');
%         tofit = medfilt1(tofit);
        tofit = tofit./max(tofit(:));

%         tofit = tofit(2:end);
        
        x = secdsback;
%         x = secdsback(2:end);

        
        fun = @(params, t)params(3).*((1+params(1)*params(2).*t).^(-params(2).^-1));
        f = fit(x, (tofit), 'exp1');
        param0 = [2, 0.01, 1];
        fitted_param = lsqcurvefit(fun, param0, x, (tofit));
        y = f.a.*exp(f.b.*x);
        yfull = f.a.*exp(f.b.*secds);
        fitted_param = real(fitted_param);
        ney = fitted_param(3).*((1+fitted_param(1)* real(fitted_param(2)) .*x).^(- real(fitted_param(2)) .^-1));
        neyfull = fitted_param(3).*((1+fitted_param(1)* real(fitted_param(2)) .*secds).^(- real(fitted_param(2)) .^-1));
        
%         RSSney = 100*sum((tofit - ney).^2)./dimension(3)
%         RSSy = 100*sum((tofit - y).^2)./dimension(3)

        tcs(1, i) = abs(f.b);
        tcs(2, i) = fitted_param(1);
        tcs(3, i) = fitted_param(2);
        
        str1 = '$$ \mathrm{Normalized} \hspace{2pt} \mathrm{data} $$';
        str2 = '$$ I(t) = A_{0}e^{kt} $$';
        str3 = '$$ I(t) = A_{0}(1 + kct)^{-\frac{1}{c}} $$';
        
        
%         hold on;
%         plot(x, tofit, 'color', col(:, i), 'LineStyle', '--', 'DisplayName', str1);

%         plot(secds, medfilt1(medfilt1(tofitfull)), 'color', col(:, i), 'LineStyle', '-.', 'DisplayName', str1);
%         plot(x, (tofitfull), 'color', col(:, i), 'LineStyle', '-.', 'DisplayName', str1);
        plot(x, (tofit), 'color', col(:, i), 'LineStyle', '-');

%         plot(x, y, 'color', col(:, i), 'LineStyle', ':', 'DisplayName', str2);
        plot(x, ney, 'color', col(:, i), 'LineStyle', '-.');

%         plot(secds2, y2, 'color', col(:, i), 'LineStyle', '-.');

    catch
        i
    end
%     plot(secds, squeeze(Peaks_all(i, 3, :))./max(squeeze(Peaks_all(i, 3, :))), 'color', col(:, i), 'LineStyle', ':', 'DisplayName', strcat('raw data', num2str(i)));
%     plot(secds, conv(squeeze(Peaks_all(i, 3, :)), fspecial('gaussian', [10,1], 5), 'same'), 'color', col(:,i), 'LineStyle', '-', 'DisplayName', strcat('raw data', num2str(i)));

    hold off;

end

title(inputname(1));
h = legend(str1, str3);
set(h,'Interpreter','latex');
set(h,'FontSize',18);
ylabel('$$ \frac{\mathrm{Mean \hspace{2pt} intensity \hspace{2pt} over \hspace{2pt} spot \hspace{2pt} area}}{\mathrm{Mean \hspace{2pt} intensity \hspace{2pt} of \hspace{2pt} background}} $$', 'Interpreter','LaTex', 'FontSize',18);
xlabel('$$ t \mathrm{(s)} $$', 'Interpreter','LaTex', 'FontSize',18);
legend('show');


figure(45)
clf();
plot(x, squeeze(Peaks_all(num_peaks+1, 3, :)));


figure(90)
clf();
for i = 1:num_peaks;
    hold on;
%     toplot = (conv(log(squeeze(Peaks_all(i, 3, :))), fspecial('gaussian', [5,1], 3), 'valid'));
    toplot = log(squeeze(Peaks_all(i, 3, :)));
%     toplot2 = log(squeeze(Peaks_all(i, 3, :)));    
    secds1 = secdsback(1:length(toplot));
%     secds2 = secdsback(1:(length(toplot2)));

    plot(secds1, toplot, 'color', col(:, i), 'LineStyle', '-', 'DisplayName', strcat('log of gaussian', num2str(i)));
%     plot(secds2, toplot2, 'color', col(:, i), 'LineStyle', '-.', 'DisplayName', strcat('log of gaussian', num2str(i)));

end
hold off;


1./tcs(1, :)

threads_start = [];
threads = [];


% 
% figure(39);
% clf();
% plot(squeeze(Peaks_all(num_peaks+1, 3, :)), 'color', 'b');




end

