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
for signal_amplitude = [0.00:delta:0.5]        
        for SpatialFreq = [10]
%                 disp(SpatialFreq)
            for eccentricity = [0]
%                 disp(eccentricity)
                grat = zeros(imsize_x,imsize_y);
%                 mask = zeros(imsize_x,imsize_y);
%                 im = zeros(imsize_x,imsize_y);
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
                        grat(ii, jj) = scale;
                    end
                end  
            end
        end

    % noise1= ((noise_amplitude*((noise - .5) * 2)) + noise_amplitude)/2;
    % grat1 = ((signal_amplitude*grat)+signal_amplitude)/2;
    %noise1 = ((noise_amplitude*((noise - .5) * 2)) + 1)/2;

    noise1 = (noise - 0.5) * 2;
    noise2 = noise1 * noise_amplitude; 
    min_noise = min(min(noise2));
    max_noise = max(max(noise2));
    fprintf('\nMin Noise %f\n', min_noise);
    fprintf('Max Noise %f\n',max_noise);

    grat1 = signal_amplitude*grat;
    min_grat = min(min(grat1));
    max_grat = max(max(grat1));
    fprintf('\nMin Grat %f\n', min_grat);
    fprintf('Max Grat %f\n',max_grat);

    %To check number of gray levels in the grating as a function of amp
%     grat2 = zeros(imsize_x,imsize_y);
    grat2 = round(((grat1/2)+0.5)*255);
    levels(1,count) = length(unique(grat2));
    levels(2,count) = signal_amplitude;
    levels(3,count) = median(median(grat2));
    levels(4,count) = max(max(grat2));
    levels(5,count) = min(min(grat2));

    im1 = noise2+grat1;
    im2 = (im1/2)+0.5;
%     figure(2)
%     imagesc(im2);colorbar ; colormap gray; set(gca,'CLim',[0,1]);
    min_image = min(min(im2));
    max_image = max(max(im2));
    fprintf('\nMin Image %f\n', min_image);
    fprintf('Max Image %f\n',max_image);
    % test = (grat1(:, 960)/2)+0.5;
    % % test = im2(:, 960);
    % figure(2);
    % plot(test);
    count = count+1;
end
plot([levels(2,:)],[levels(4,:)],'d-'); hold on
plot([levels(2,:)],[levels(5,:)],'d-'); hold on
xlabel('Signal Amp')
ylabel('Gray levels')  

if any(any(diff(levels(4:5,:)') == 0))
    title('bad step')
end
% plot([levels{2,:}],[levels{1,:}])
% xlabel('Signal Amp')
% ylabel('Gray levels')