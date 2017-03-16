function [ posRefine ] = Gauss_fit_to_peak( I, pos )


% This function refines 2D peaks using a 5 parameter Gaussian:
%  I = I0 * exp( - ((x-x0).^2 + (y-y0).^2) / (2*sigma^2) ) + k
% I0 is the peak intensity, x0 and y0 are the refined peak positions in
% pixels, sigma is the 2D Gaussian standard deviation, k is the background.

% Inputs:
% I         - 2D image array.  Should probably be doubles or singles.
% pos       - N x 2 array containing peak positions to be refined.

% Outputs:
% posRefine - N x 5 array containing refined positions.
% [x_fit y_fit I0_fit sigma_fit k_fit]

% % Refine peak positions with nonlinear least squares 2D Gaussian fit.
% % Fitting variables

% STEM fits
NsubPxIter = 16;        % Sub pixel fitting iterations.
rCut = 8;               % Cutting radius (to speed up fits).
rFit = 4;               % fitting radius in pixels.
sigma0 = 2;             % initial guess for 2D Gaussian sigma .
xyStep = 1;           % maximum change in x/y per iteration.
sStep = 1/2;            % maximum change in sigma per iteration.
f_plot = 1;             % To plot as you go, set to "true" or 1.
sigmaRange = [0.5 4];    % minimum and maximum allowed sigma values.

flag_quiet = 0;   % Set to 1 to suppress console completion messages.

fun1 = @(c,x) c(1)*exp((-1/2/c(2)^2) ...
    *((x(:,1)-c(3)).^2+(x(:,2)-c(4)).^2)) + c(5);
options = optimset('TolX',1e-3,'TolFun',1e-3,...
    'MaxFunEvals',10^6);

% coordinates
N = double(size(I));
[ya,xa] = meshgrid(1:N(2),1:N(1));
v = -rCut:rCut;
r2 = rFit^2;

Np = size(pos,1);
posRefine = zeros(Np,5);

%Stepping over each peak in the list:
for a0 = 1:Np
    % Cut out segment
    x = pos(a0,1);
    y = pos(a0,2);
    
    if x > 1 + rCut && x < N(1)-rCut ...
            && y > 1 + rCut && y < N(2)-rCut
        
        xv = v + round(x);
        yv = v + round(y);
        xv(xv<1) = [];
        yv(yv<1) = [];
        xv(xv>N(1)) = [];
        yv(yv>N(2)) = [];
        
        
        xCut = xa(xv,yv);
        yCut = ya(xv,yv);
        ICut = double(I(xv,yv));
        
        % Estimate parameters
        k = min(ICut(:));
        I0 = max(ICut(:)) - k;
        s = sigma0;
        
        for a1 = 1:NsubPxIter
            % Subset
            sub = (xCut-x).^2 + (yCut-y).^2 < r2;
            xfit = xCut(sub);
            yfit = yCut(sub);
            zfit = ICut(sub);
            
            % Initial guesses and bounds
            c0 = [I0 s x y k];
            dI = 0.2*I0;
            lb = [max(I0-dI,0) max(s-sStep,sigmaRange(1)) x-xyStep y-xyStep k-dI];
            ub = [I0+dI min(s+sStep,sigmaRange(2)) x+xyStep y+xyStep k+dI];

            try
                %Perform fitting
                [~,cc] = evalc('lsqcurvefit(fun1,c0,[xfit yfit],zfit,lb,ub,options)');
            catch ME
                continue
            end
            
            I0 = cc(1);%*(1-damp) + I0*damp;
            s = cc(2);
            x = cc(3);%*(1-damp) + x*damp;
            y = cc(4);%*(1-damp) + y*damp;
            k = cc(5);
            
            %             disp(num2str([I0 x y]))
            %             Itest = I0*exp((-1/2/sMean^2)*((xCut-x).^2+(yCut-y).^2));
            if f_plot == 1
                Itest = I0*exp((-1/2/s^2)*((xCut-x).^2+(yCut-y).^2))+k;
                figure(1)
                clf
                cR = [min(ICut(sub)) max(ICut(sub))];
                imagesc([ICut.*sub Itest.*sub])
                axis equal off
                colormap(hot(256))
                caxis(cR)
            end
        end
        posRefine(a0,1:5) = [x y I0 s k];
        
        if flag_quiet ~= 1
            comp = a0 / Np;
            if mod(a0,100) == 0
                disp(['Fitting is ' sprintf('%.02f',comp*100) '% complete'])
            end
        end
        
    else
        posRefine(a0,1:5) = [x y 0 0 0];
    end
    %progressbar(comp,2);
    
end


end

