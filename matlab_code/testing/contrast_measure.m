noise_path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Matlab_Code\main\noise';
directory = dir(noise_path);

imsize = 1920;
% pixelAngle = 0.8; % arcmin per pixel
pixelAngle = 1.06; %Vis dynamics

width = 1;
gaussSD = .2;
phase = 0;

noise_amplitude = .7;
% amplitude = .3;
M_PI = pi;

% eccentricity = 2;
% SpatialFreq = [2 10];
contrast = struct;
% for i=3:4
for i=3:length(directory)
    filePath = fullfile(noise_path,directory(i).name);
    noise = imread(filePath); 
    count=1;
%     imshow(noise)
%     pause
%     close all
    for amplitude = 0.01:0.01:0.3        
        for SpatialFreq = [2 10]
%                 disp(SpatialFreq)
            for eccentricity = [0 4 8]
%                 disp(eccentricity)
                grat = zeros(imsize);
                mask = zeros(imsize);
                im = zeros(imsize);
                for ii = 1:imsize
                    for jj = 1:imsize
                        r = sqrt ((ii - imsize/2)^2 + (jj - imsize/2)^2) * pixelAngle / 60;
                        mask(ii,jj) = false;
                        if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2))
                            scale = amplitude * sin(2 * M_PI * SpatialFreq * r + phase);
                            mask(ii,jj) = true;
                        elseif (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD)
                            scale = amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) *...
                                (cos(M_PI / gaussSD * (r - eccentricity - width / 2)) + 1) / 2;
                        elseif (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD)
                            scale = amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) *...
                                (cos(M_PI / gaussSD * (r - eccentricity + width / 2)) + 1) / 2;
                        else
                            scale = 0;
                        end                       
                        
                        grat(ii, jj) = scale;
                    end
                end  
                
                noise = (noise - .5) * noise_amplitude + .5;
%                 noise = double(noise);
%                 grat = uint8(grat);
%                 grat = repmat(grat,[1 1 3]);
                noise = double(noise(:,:,1));
                im = noise + grat;
                im(~mask) = nan;
%                 imagesc(im);
                michelson_con = (max(max(im)) - min(min(im)))/( min(min(im)) + max(max(im)));
%                 pause
                contrast(i-2).amp(count) = amplitude;
                contrast(i-2).freq(count) = SpatialFreq;
                contrast(i-2).ecc(count) = eccentricity;
                contrast(i-2).con(count) = michelson_con;
                count = count+1;
            end         

        end   
    end
end
low_amp = [];
for i=1:10
    low_amp = [low_amp;contrast(i).con(1:180)];
end