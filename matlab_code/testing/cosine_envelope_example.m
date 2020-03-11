% 1 dimensional example of cosine-enveloped radial wave image
% Janis Intoy, April 17, 2019

% parameters all in units of arcmin
SpatialFreq = 2 / 60;
eccentricity = 4 * 60;
width = 1.5 * 60;
cos_halfperiod = 1 * 60; % this is called gaussSD in the shader
phase = 0;

M_PI = pi;

rspace = linspace(0, 8*60, 500);

img = nan(size(rspace));
env = nan(size(rspace));

gaussSD = cos_halfperiod; % this is its name in the shader (remnant of gaussian envelope)
for ri = 1:length(rspace)
    r = rspace(ri);
    
    %%%% this block is almost exactly copy-pasted from the shader %%%%%
    if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2))
        signal =  sin(2 * M_PI * SpatialFreq * r + phase);
    elseif (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD)
        signal =  sin(2 * M_PI * SpatialFreq * r  + phase)* ...
        (cos(M_PI / gaussSD * (r - eccentricity - width / 2))+ 1) / 2;
        
    elseif (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD)
        signal = sin(2 * M_PI * SpatialFreq * r  + phase)* ...
        (cos(M_PI / gaussSD * (r - eccentricity + width / 2))+ 1) / 2;
        
    else
        signal = 0;
    end
    %%%% this block is almost exactly copy-pasted from the shader %%%%%
    
    % this part is just the raised cosine envelope part of signal
    if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2))
        env(ri) = 1;
    elseif (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD)
        
        env(ri) = (cos(M_PI / gaussSD * (r - eccentricity - width / 2))+ 1) / 2;
        
    elseif (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD)
        env(ri) = (cos(M_PI / gaussSD * (r - eccentricity + width / 2))+ 1) / 2;
        
    else
        env(ri) = 0;
    end
    
    img(ri) = signal;
end

%% plot 1d image
figure(1); clf; hold on;
% 'raised' portion of envelope
patch([eccentricity-width/2, eccentricity-width/2, eccentricity+width/2, eccentricity+width/2],...
    [-1, 2, 2, -1],...
    'g', 'EdgeColor', 'none', 'FaceAlpha', .3);
% cosine envelope on either side
patch([eccentricity-width/2, eccentricity-width/2, eccentricity-width/2 - cos_halfperiod, eccentricity-width/2-cos_halfperiod],...
    [-1, 2, 2, -1],...
    'r', 'EdgeColor', 'none', 'FaceAlpha', .3);
patch([eccentricity+width/2, eccentricity+width/2, eccentricity+width/2 + cos_halfperiod, eccentricity+width/2+cos_halfperiod],...
    [-1, 2, 2, -1],...
    'r', 'EdgeColor', 'none', 'FaceAlpha', .3);

hp1 = plot(rspace, img, 'b', 'linewidth', 2);
hp2 = plot(rspace, env, 'k', 'linewidth', 2);

xlabel('eccentricity (arcmin)');
set(gca, 'FontSize', 14);

xlim([100, 400]);
ylim([-1, 2]);

text(eccentricity,...
    1.7, sprintf('annulus width\n= %i arcmin', width),...
    'HorizontalAlignment', 'center', 'FontSize', 12);

% cosine envelope
text(mean([eccentricity-width/2 - cos_halfperiod, eccentricity-width/2]),...
    1.7, sprintf('cos half period\n= %i arcmin', cos_halfperiod),...
    'HorizontalAlignment', 'center', 'FontSize', 12);
text(mean([eccentricity+width/2 + cos_halfperiod, eccentricity+width/2]),...
    1.7, sprintf('cos half period\n= %i arcmin', cos_halfperiod),...
    'HorizontalAlignment', 'center', 'FontSize', 12);
errorbar([eccentricity-width/2 - cos_halfperiod, eccentricity-width/2],...
    [1.5, 1.5], [.1, .1], 'k', 'linewidth', 2, 'capsize', 0);
errorbar([eccentricity+width/2 + cos_halfperiod, eccentricity+width/2],...
    [1.5, 1.5], [.1, .1], 'k', 'linewidth', 2, 'capsize', 0);

legend([hp1, hp2], {'enveloped signal', 'raised cosine envelope'}, 'Location', 'southeast');
