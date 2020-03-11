imsize_x = 1080;
imsize_y = 1920;
pixelAngle = 1.06; % arcmin per pixel, was 0.8 initially
width = 1;
gaussSD = 0.2;
phase = 0;
M_PI = pi;
noise_amplitude = 0.6;
noise = nnoise_nk(imsize_x,imsize_y, pixelAngle);
noise = ((noise / max(abs(noise(:)))) / 2) + 0.5;%range(0,1)

count = 1;
delta = 0.01;%.008 looks good
% signal_amplitude = 0.3;
eccentricity = 4;
SpatialFreq = 10;
grat = zeros(imsize_x,imsize_y,3);

% for signal_amplitude = [0.00:delta:0.5]
for signal_amplitude = 0.3
%     for channel = 1:1:3
        for ii = 1:imsize_x
            for jj = 1:imsize_y
                r = sqrt ((ii - imsize_x/2)^2 + (jj - imsize_y/2)^2) * pixelAngle / 60;
                %                         mask(ii,jj) = false;
                if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2))
                    scale = sin(2 * M_PI * SpatialFreq * r + phase);
                    %                             mask(ii,jj) = true;
                elseif (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD)
                    scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                        (cos(M_PI / gaussSD * (r - eccentricity - width / 2))+ 1) / 2;
                elseif (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD)
                    scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                        (cos(M_PI / gaussSD * (r - eccentricity + width / 2))+ 1) / 2;
                else
                    scale = 0;
                end
                grat(ii, jj, 1) = scale;
                grat(ii, jj, 2) = scale;
                grat(ii, jj, 3) = scale;
            end
        end
%     end
    noise1 = (noise - 0.5) * 2;
    noise2 = noise1 * noise_amplitude;
%     min_noise = min(min(noise2));
%     max_noise = max(max(noise2));
%     fprintf('\nMin Noise %f\n', min_noise);
%     fprintf('Max Noise %f\n',max_noise);
    
    grat1 = signal_amplitude*grat;
%     min_grat = min(min(grat1));
%     max_grat = max(max(grat1));
%     fprintf('\nMin Grat %f\n', min_grat);
%     fprintf('Max Grat %f\n',max_grat);
    
    %To check number of gray levels in the grating as a function of amp
    %     grat2 = zeros(imsize_x,imsize_y);
    grat2 = round(((grat1/2)+0.5)*255);
    levels(1,count) = length(unique(grat2));
    levels(2,count) = signal_amplitude;
%     levels(3,count) = median(median(grat2));
%     levels(4,count) = max(max(grat2));
%     levels(5,count) = min(min(grat2));
    
    im1 = noise2+grat1;
    im2 = (im1/2)+0.5;
    
    figure(1)
    imagesc(im2);
    colorbar ; 
    colormap gray; 
    set(gca,'CLim',[0,1]);
    
    min_image = min(min(im2));
    max_image = max(max(im2));
    fprintf('\nOriginal\n');
    fprintf('\nMin Image %f\n', min_image);
    fprintf('Max Image %f\n',max_image);
    % test = (grat1(:, 960)/2)+0.5;
    % % test = im2(:, 960);
    % figure(2);
    % plot(test);
    count = count+1;
    
end

%% Bit-stealig + Dithering
grat_dith = zeros(imsize_x,imsize_y,3);
LUT = [0,0,0, 0,0,1, 1,0,0,	1,0,1, 0,1,0, 0,1,1, 1,1,0,	1,1,1];
im3 = im2*255;
for ii = 1:imsize_x
    for jj = 1:imsize_y
        grayImg = ((im3(ii,jj,1) * 2)+ (im3(ii,jj,2) * 4)+(im3(ii,jj,3) * 1))/7; %//R:G:B -> 2:4:1
%         grayImg = im3(ii,jj);
        intPart = floor(grayImg);
        step = (grayImg - intPart)*7;
        stepInt = floor(step);
        stepDec = step - stepInt;
        rnd = rand();
        
        if (rnd <  stepDec )
            stepRnd = stepInt + 1;           
        elseif (rnd >= stepDec )
            stepRnd = stepInt;            
        end
        
        grat_dith(ii, jj, 1) = intPart + LUT(3 * stepRnd + 1);
        grat_dith(ii, jj, 2) = intPart + LUT(3 * stepRnd + 2);
        grat_dith(ii, jj, 3) = intPart + LUT(3 * stepRnd + 3);
        
    end
end
grat_dith_scaled = grat_dith/255;
figure(2)
imagesc(grat_dith_scaled);
colorbar ;
colormap gray;
set(gca,'CLim',[0,1]);

fprintf('\nDithered\n');
min_image = min(min(grat_dith_scaled));
max_image = max(max(grat_dith_scaled));
fprintf('\nMin Image %f\n', min_image);
fprintf('Max Image %f\n',max_image);

