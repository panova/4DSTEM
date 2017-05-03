function [ Grain_map, full_grain_array  ] = Grain_separation_flooding( FFA, realDimX, realDimY, mycmap )
%GRAIN_SEPARATION_FLOODING Summary of this function goes here
%   Detailed explanation goes here

dimension = size(FFA);
padd = 2;

FFA = reshape(FFA, [dimension(1), realDimX, realDimY]);
FFA = padarray(FFA, [0, padd, padd], 'both');
size(FFA);

Grain_map = zeros([realDimX, realDimY]);
Grain_map = padarray(Grain_map,[padd, padd], -1, 'both');
Grain_map_size = size(Grain_map);

Candidates_map = zeros([realDimX, realDimY]);
Candidates_map = padarray(Candidates_map,[padd, padd], -1, 'both');

del_angle_direct = 3;
del_angle_diagonal = 3;

full_grain_array = zeros([Grain_map_size(1), Grain_map_size(2), 124, 124]);

% for C = 2:realDimY-1
%This is the plotted X on imagesc
for globalC = 3:126; %55
    globalC
    %     for R = 2:realDimX-1
    %This is the plotted Y on imagesc
    for globalR = 3:126 %82
%         pause(0.05);
        fig = 5;
        Grain_map = zeros([realDimX, realDimY]);
        Grain_map = padarray(Grain_map,[padd, padd], -1, 'both');
        
        Candidates_map = zeros([realDimX, realDimY]);
        Candidates_map = padarray(Candidates_map,[padd, padd], -1, 'both');

        tracker = zeros([realDimX, realDimY]);
        tracker = padarray(tracker,[padd, padd], -1, 'both');
        tracker(globalR, globalC) = 1;

        peaks = FFA(:, globalR, globalC);
        peak = peaks(4);
        if peak == -1;
            'There is no diffraction here!'
            continue
        end
        Grain_map(globalR, globalC) = peak;
        Candidates_map(globalR, globalC) = 1;

        
        rows = globalR;
        cols = globalC;
        
        cnt = 1;
        
        length(Candidates_map(Candidates_map == 1));
        
        while length(Candidates_map(Candidates_map == 1)) > 0 && cnt < 300
        
%         cnt
%         while cnt < 5;

        for candidateID = 1:length(rows)
            C = cols(candidateID);
            R = rows(candidateID);

                peak = Grain_map(R, C);
                %Direct neighbors:
                %LEFT
                if Grain_map(R, C-1) == 0
                    pleft = FFA(:, R, C-1);
                    if length(pleft(abs(pleft -  peak) < del_angle_direct)) >= 1
                        Grain_map(R, C-1) = mean(pleft(abs(pleft -  peak) < del_angle_direct));
                        Candidates_map(R, C-1) = 1;
                    else
                        Grain_map(R, C-1) = -1;
                        Candidates_map(R, C-1) = -1;
                    end
                end
                
                %RIGHT
                if Grain_map(R, C+1) == 0
                    pright = FFA(:, R, C+1);
                    if length(pright(abs(pright -  peak) < del_angle_direct)) >= 1
                        Grain_map(R, C+1) = mean(pright(abs(pright -  peak) < del_angle_direct));
                        Candidates_map(R, C+1) = 1;
                    else
                        Grain_map(R, C+1) = -1;
                        Candidates_map(R, C+1) = -1;
                    end
                end
                
                %UP
                if Grain_map(R-1, C) == 0
                    pup = FFA(:, R-1, C);
                    if length(pup(abs(pup -  peak) < del_angle_direct)) >= 1
                        Grain_map(R-1, C) = mean(pup(abs(pup -  peak) < del_angle_direct));
                        Candidates_map(R-1, C) = 1;
                    else
                        Grain_map(R-1, C) = -1;
                        Candidates_map(R-1, C) = -1;
                    end
                end
                
                %DOWN
                if Grain_map(R+1, C) == 0
                    pdown = FFA(:, R+1, C);
                    if length(pdown(abs(pdown -  peak) < del_angle_direct)) >= 1
                        Grain_map(R+1, C) = mean(pdown(abs(pdown -  peak) < del_angle_direct));
                        Candidates_map(R+1, C) = 1;
                    else
                        Grain_map(R+1, C) = -1;
                        Candidates_map(R+1, C) = -1;
                    end
                end
                
                
                
                %Diagonal neighbors:
                
                %LEFT UP
                if Grain_map(R-1, C-1) == 0
                    pleftup = FFA(:, R-1, C-1);
                    if length(pleftup(abs(pleftup -  peak) < del_angle_diagonal)) >= 1
                        Grain_map(R-1, C-1) = mean(pleftup(abs(pleftup -  peak) < del_angle_diagonal));
                        Candidates_map(R-1, C-1) = 1;
                    else
                        Grain_map(R-1, C-1) = -1;
                        Candidates_map(R-1, C-1) = -1;
                    end
                end
                
                %RIGHT UP
                if Grain_map(R-1, C+1) == 0
                    prightup = FFA(:, R-1, C+1);
                    if length(prightup(abs(prightup -  peak) < del_angle_diagonal)) >= 1
                        Grain_map(R-1, C+1) = mean(prightup(abs(prightup -  peak) < del_angle_diagonal));
                        Candidates_map(R-1, C+1) = 1;
                    else
                        Grain_map(R-1, C+1) = -1;
                        Candidates_map(R-1, C+1) = -1;
                    end
                end
                
                %LEFT DOWN
                if Grain_map(R+1, C-1) == 0
                    pleftdown = FFA(:, R+1, C-1);
                    if length(pleftdown(abs(pleftdown -  peak) < del_angle_diagonal)) >= 1
                        Grain_map(R+1, C-1) = mean(pleftdown(abs(pleftdown -  peak) < del_angle_diagonal));
                        Candidates_map(R+1, C-1) = 1;
                    else
                        Grain_map(R+1, C-1) = -1;
                        Candidates_map(R+1, C-1) = -1;
                    end
                end
                
                %RIGHT DOWN
                if Grain_map(R+1, C+1) == 0
                    prightdown = FFA(:, R+1, C+1);
                    if length(prightdown(abs(prightdown -  peak) < del_angle_diagonal)) >= 1
                        Grain_map(R+1, C+1) = mean(prightdown(abs(prightdown -  peak) < del_angle_diagonal));
                        Candidates_map(R+1, C+1) = 1;
                    else
                        Grain_map(R+1, C+1) = -1;
                        Candidates_map(R+1, C+1) = -1;
                    end
                end
                
                Candidates_map(R, C) = -1;
                
%                 alpha = (~(Grain_map <= 0))*1;
%                 figure(97);
%                 clf();
%                 imagesc(Candidates_map);
%                 caxis([min(Candidates_map(:)), max(Candidates_map(:))]);
%                 title('New Candidates');
%                 axis equal off
%                 colormap(violetFire(256))
%                 drawnow;
%                 
%                 
%                 figure(98);
%                 clf();
%                 imagesc(Grain_map, 'AlphaData', alpha);
% %                 caxis([min(Grain_map(Grain_map>0)), max(Grain_map(Grain_map>0))]);
%                 caxis([0 180])
%                 title('Grain Map');
%                 colormap(mycmap);
%                 axis equal off
%                 drawnow;

            
        end
        
        length(Candidates_map(Candidates_map == 1));
        [rows, cols, ~] = find(Candidates_map == 1);
        cnt = cnt +1;
        
        
        end
%         
%         alpha = (~(Grain_map <= 0))*1;
%         figure(fig+1);
%         clf();
%         imagesc(Grain_map, 'AlphaData', alpha);
%         %caxis([min(Grain_map(Grain_map>0)), max(Grain_map(Grain_map>0))]);
%         caxis([0 180])
%         title('Grain Map');
%         colormap(mycmap);
%         axis equal off
%         drawnow;
%         
%         figure(3469)
%         clf();
%         imagesc(tracker);
%         axis equal off;
        
        full_grain_array(:, :, R, C) = Grain_map;
    end
    
end

% 
% alpha = (~(Grain_map <= 0))*1;
% figure(fig);
% clf();
% imagesc(Candidates_map);
% caxis([min(Candidates_map(:)), max(Candidates_map(:))]);
% title('New Candidates');
% axis equal off
% colormap(violetFire(256))
% drawnow;

% 
% figure(fig+1);
% clf();
% imagesc(Grain_map, 'AlphaData', alpha);
% %caxis([min(Grain_map(Grain_map>0)), max(Grain_map(Grain_map>0))]);
% caxis([0 180])
% title('Grain Map');
% colormap(mycmap);
% axis equal off
% drawnow;




end

