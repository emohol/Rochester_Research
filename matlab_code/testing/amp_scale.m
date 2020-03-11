% imsize = 1920;
imsize_x = 1080;
imsize_y = 1920;
pixelAngle = 1.06; % arcmin per pixel, was 0.8 initially
eccentricity = 4;
SpatialFreq = 2;
M_PI = pi;
width = 1;
gaussSD = .2;
phase = 0;
% noise_amplitude = .5;
% signal_amplitude = .45;

% raw_noise=rand(imsize);
% noise1 = noise_amplitude*((raw_noise*2)-1);

noise_orig = nnoise_nk(imsize_x,imsize_y, pixelAngle);

% noise_scaled = noise_orig / max(abs(noise_orig(:))); % range(-1,1)
noise_scaled = noise_orig / max(abs(noise_orig(:))) / 2 + 1/2;%range(0,1)
% fprintf('Min Pink Noise %f\n',min(min(noise_scaled)))
% fprintf('Max Pink Noise %f\n\n',max(max(noise_scaled)))
% noise_scaled = noise_amplitude*noise_scaled;
% fprintf('Min Scaled Pink Noise %f\n',min(min(noise_scaled)))
% fprintf('Max Scaled Pink Noise %f\n\n',max(max(noise_scaled)))
% noise1 = (noise_orig*2)-1;
% noise = noise_amplitude*noise_orig;
A = [];
count = 1;
delta_sig = .01;
% noise_path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Matlab_Code\main\noise_nk(pos)';
% directory = dir([noise_path filesep '*.bmp']);
% for i=1:length(directory)
%     filePath = fullfile(noise_path,directory(i).name);
%     noise = double(noise_orig(:,:,1));
    for noise_amplitude = [0.05:.05:.95]   
        noise1 = (noise_scaled - 0.5) * 2;
        noise = noise1 * noise_amplitude; 
        fprintf('Min Scaled Pink Noise %f\n',min(min(noise)))
        fprintf('Max Scaled Pink Noise %f\n\n',max(max(noise)))
        A(1,count)=noise_amplitude;
        for signal_amplitude = [0:delta_sig:.99]   
            grat = zeros(imsize_x,imsize_y);
            for ii = 1:imsize_x
                for jj = 1:imsize_y
                    r = sqrt ((ii - imsize_x/2)^2 + (jj - imsize_y/2)^2) * pixelAngle / 60;
                    if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2))
                        scale = sin(2 * M_PI * SpatialFreq * r + phase);
                    elseif (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD)
                        scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                            (cos(M_PI / gaussSD * (r - eccentricity - width / 2)) + 1) / 2;
                    elseif (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD)
                        scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                            (cos(M_PI / gaussSD * (r - eccentricity + width / 2)) + 1) / 2;
                    else
                        scale = 0;
                    end
                    grat(ii, jj) = scale;
                end
            end
            grat1 = signal_amplitude*grat;
            fprintf('Min Grating %f\n', min(min(grat1)));
            fprintf('Max Grating %f\n\n',max(max(grat1)))
    %         img = zeros(imsize_x,imsize_y);
            img = noise + grat1; %both range (-1,+1), this is when you add the two **IMPORTANT
            image = (img/2)+0.5; %range (0,1). Note that you divide first when range (-1,1) then add
            fprintf('Min Image %f\n',min(min(image)))
            fprintf('Max Image %f\n\n',max(max(image)))
            if any(any(image<0))
                warning('negative values in image, re-range [0 1], amp = %0.2f', signal_amplitude)
                A(2,count)=signal_amplitude-delta_sig;
                break;
            elseif any(any(image>1))
                warning('extreme values in image, re-range [0 1],  amp = %0.2f', signal_amplitude)
                A(2,count)=signal_amplitude-delta_sig;
                break;
            end
            % imshow(image,[]);
            image_255 = round(255 * image);%[0,255]
        %     clear
        end
        count = count+1;
    end
% end