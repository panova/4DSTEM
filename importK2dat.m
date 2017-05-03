function [stack3D] = importK2dat()
tic
% Import K2 data cube (3D) that has been converted to binary uint16.
% If virtual detector is provided, create hologram image from det.

% flag_remove_xrays = 0;  % xray threshold, set to zero to ignore.
fname = cell(3,1);

%A1; lattice_rad = 63, center at [125 120]
% fname = ['C:\Users\NCEM\Documents\MATLAB\073015_4DSTEM_P3HT\P3HTNW_AN_50x50_40000_CL=245_a=p51_spot11\' ...
%     'P3HTNW_AN_50x50_40000_CL=245_a=p51_spot11_CROPPED.dat']; 
% xySize = [253 266];

%A2 lattice_rad = 63, same as A1
% fname = ['C:\Users\NCEM\Documents\MATLAB\073015_4DSTEM_P3HT\P3HTNW_ASCAST_50x50_40000_CL=245_a=p51_spot11\' ...
%     'P3HTNW_ASCAST_50x50_40000_CL=245_a=p51_spot11_CROPPED.dat']; 
% xySize = [253 265];

% %A3 lattice_rad = 130
% fname = ['C:\Users\NCEM\Documents\MATLAB\073015_4DSTEM_P3HT\P3HTAN_50x50_40000_CL=480_a=p51_spot11_1\' ...
%     'P3HTAN_50x50_40000_CL=480_a=p51_spot11_1.dat']; 
% xySize = [512 512];

%D2
% fname = ['C:\Users\NCEM\Documents\MATLAB\100905_4DSTEM_EDS_P3HT_PS\5_5_annealed1h_200kV\' ...
%     'P3HT_PS_5_5_annealed 1h_128x128_20000.dat']; 
% xySize = [512 512];

%D1
% fname = ['C:\Users\NCEM\Documents\MATLAB\100905_4DSTEM_EDS_P3HT_PS\5_5_annealed1h_200kV\' ...
%     'P3HT_PS_5_5_annealed 1h_64x64_40000.dat']; 
% xySize = [512 512];

% %C1
% fname = ['C:\Users\NCEM\Documents\MATLAB\100615_4DSTEM_EDS_P3HT_PS\' ...
%     '5_5 P3HT_PS_20000_50x50_CROPPED.dat']; 
% xySize = [232 239];
% 

%E1
% fname = ['C:\Users\NCEM\Documents\MATLAB\11102015_P3HTPSS_50_50_TITAN\' ...
%     '11102015_P3HTPSS_50_50_CL245mm_spot11_64x64_40nm.dat']; 
% xySize = [512 512];

%PEO
% fname = ['C:\Users\NCEM\Documents\MATLAB\PEO\' ...
%     'PEO10K_THICK_CL756_ID3_128x128_40nmstep_bin8_0p1s.dat'];


%PE
% fname = ['C:\Matlab\Ouliana Panova\NanoDiffraction\PE\Winey TITAN PE 2016\' ...
%     'KCB HDPE 20160120 Stack 22 33ms spot 10 ss70nm 60x60 alpha=p51mrad_cropped.dat']; 
% xySize = [512 512];

%PE WINEY KAREN K2 0-02-2016, #14
% fname = ['C:\Data\Ouliana\4DSTEM\PE WINEY KAREN K2 0-02-2016\' ...
%     'r45AA_Xtal1_14_stack.dat']; 
% xySize = [506 410];

%PE WINEY KAREN K2 0-02-2016, #3
% fname = ['C:\Data\Ouliana\4DSTEM\PE WINEY KAREN K2 0-02-2016\' ...
%     'r45AA_Xtal1_3_stack.dat']; 
% xySize = [485 409];

%PE WINEY KAREN K2 0-02-2016, #18
% fname = ['C:\Data\Ouliana\4DSTEM\PE WINEY KAREN K2 0-02-2016\' ...
%     'r45AA_Xtal1_18_stack.dat']; 
% xySize = [489 410];

% % P3HT:PS 50:50 CRYO K2 August 28th 2016 MM12
% fname = ['E:\Data\Ouliana\12\' ...
%     'p3ht12_Processed_stack.dat']; 
% 
% xySize = [488 488];

% %20160928 HDPE Sample ID 2 LN2 Karen 
% fname = ['E:\Data\Ouliana\20160928 HDPE Sample ID 2\' ...
%     'Stack 5 64x64 80ss 0p48 CL=480 bin4 33ms 300kV spot 11_180C.dat']; 
% 
% xySize = [512 512];

%PEO 4K thick TITAN 2 nm
% fname = ['E:\Data\Ouliana\08142016TITAN_PEO_4K_thick\' ...
%     'Stack 5 64x64 2nm step 33ms CL=480 alpha=p48 bin4 300kV.dat']; 
% 
% xySize = [512 512];

%PEO 4K thick TITAN 5 nm
% fname = ['E:\Data\Ouliana\08142016TITAN_PEO_4K_thick\' ...
%     'Stack 4 64x64 5nm step 100ms CL=480 alpha=p48 bin4 300kV.dat'];

%PEO 4K thick TITAN 5 nm
% fname = ['E:\Data\Ouliana\08142016TITAN_PEO_4K_thick\' ...
%     'Stack 3 60x60 ss=5nm 33ms CL=480 alpha=0p48 bin4 300kV great.dat'];
% 
% xySize = [512 512];

%=================================================================================================
% TAKACS TITAN CRYO no DIO sample B10 [December 8 2016] "TC7" step=10nm
% fname = ['E:\Data\Ouliana\TACAKS_Dec08_2016\' ...
%     '7TakacsB10_10082016_cryoCL380mmp033sSS10nm128x128.dat']; 

% TAKACS TITAN CRYO no DIO sample B10 [December 8 2016] "TC4" step=20nm
% fname = ['E:\Data\Ouliana\TACAKS_Dec08_2016\' ...
%     '4TakacsB10_10082016_cryoCL600mmp033sSS20nm128x128.dat']; 

% % TAKACS TITAN CRYO DIO sample C10 [January 17] "TCDIO9"
% fname = ['E:\Data\Ouliana\TACAKS_Jan17_2017_C10_andB10laterCRYO\20170117 Takas_cryo_C10\' ...
%     '9TACAKS_CL380mm_300kV_ss10nm_128x128.dat']; 

% TAKACS TITAN CRYO DIO sample C10 [January 17] "TCDIO4"
% fname = ['E:\Data\Ouliana\TACAKS_Jan17_2017_C10_andB10laterCRYO\20170117 Takas_cryo_C10\' ...
%     '4TACAKS_CL380mm_300kV_ss10nm_128x128.dat']; 

% TAKACS TITAN CRYO DIO sample B10 [JANUARY 17 2016] "TC22" step=10nm
fname = ['E:\Data\Ouliana\TACAKS_Jan17_2017_C10_andB10laterCRYO\20170117 Takas_cryo_C10\' ...
    'B102TACAKSB10_CL380mm_300kV_ss5nm_128x128.dat']; 
% 
% 
xySize = [512 512];

%=================================================================================================
% JOSE Rodrigues datasets 

% fname = ['E:\Data\Ouliana\XiPeptiodSheetsMarch7th2017\20170307 Xi\' ...
%     '15_pep99sheet_50x50_40nm=ss cl=600 alpha=0p13 300kV 10umC2 spot2 bin4 33ms.dat']; 
% 
% xySize = [512 512];

%=================================================================================================

% format = 'uint16=>uint16';
format = 'int16=>int16';

% xySize = [1920 1792];
%Dimensions of one slice within the stack
% For the K2 camera:
% xySize = [480 448];
% For the TITANX camera (Orius): 
% xySize = [512 512];
% For the LIBRA camera:
% xySize = [256 256];

% Get file size, other info
s = dir(fname);


% Calculate number of images within the stack, assume 2 bytes per pixel.
Nimages = s.bytes / prod(xySize) / 2;


%Nimages = 16512;

% Initialize
stack3D = zeros(xySize(1),xySize(2), Nimages,'int16');

% Import data
fid = fopen(fname,'r');
for a0 = 1:Nimages
    try
        %I = fread(fid,xySize,format);
%         SI  = Shift_upsampled( I, shift(a0, :) );    
%         stack3D(:,:,a0) = SI;
        stack3D(:,:,a0) = fread(fid,xySize,format);
    catch
        fclose(fid);
        error(['Warning, file: ' fname ' does not exist, or size is wrong.'])
    end
    
    
%             comp = (a0 / Nimages ...
%                 + a1 - 1) / Nfile;
%             progressbar(comp,2);
end

fclose(fid);
% if comp < 1
%     progressbar(1,2);
% end

Imean = mean(stack3D,3);

figure(1)
clf
imagesc(Imean)
axis equal off
colormap(hot(256))
%set(gca,'position',[0 0 1 1])


toc
end

% 
% function [I] = removeXrays(I,flag_remove_xrays)
% 
% 
% % Im = (circshift()) / 4;
% Im = medfilt2(I,[3 3]);
% % Im = conv2(I,uint16([1 1 1;1 0 1;1 1 1]),'same');
% sub = (I - Im) > flag_remove_xrays;
% sub([1 end],:) = false;
% sub(:,[1 end]) = false;
% I(sub) = Im(sub);
% 
% 
% end