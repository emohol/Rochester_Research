% function [michelson_con] = Untitled3(amplitude,eccentricity,SpatialFreq)

% switch nargin
%     case 1
%         SpatialFreq = 2;
%         if length(Inp)==1
%             I_0= Inp; %arbitrary fixed input
%         else
%             I = Inp;
%         end
%         V_reset = -0.070;% -70mV
%     case 2
%         I_0 = Inp;
%         V_reset = V_0;
%     otherwise
%         Inp = 0;
%         I_0 = 0.02;
%         V_reset = -0.070;% -70mV
% end

% draw a radial wave

SpatialFreq = 1;

imsize = 1920;
% pixelAngle = 0.8; % arcmin per pixel
pixelAngle = 1.06; %Vis dynamics

width = 1;
gaussSD = .2;
phase = 0;

noise_amplitude = .7;
amplitude = .3;
M_PI = pi;

eccentricity = 2;

%%
grat = zeros(imsize);
mask = zeros(imsize);
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

% figure(1); clf;
% imagesc(grat);
% axis image;

% %% power spectrum of radial wave
% spec = fftshift(fft2(grat));
% freq = 60 / pixelAngle * ((-imsize/2 + 1):imsize/2)/imsize;
% 
% figure(2); clf;
% imagesc(freq, freq, pow2db(abs(spec).^2));
% ca = caxis;
% title('power in grating');
%% power spectrum of noise
noise = nnoise(imsize, pixelAngle);
noise = noise / max(abs(noise(:))) / 2 + 1/2;

noise = (noise - .5) * noise_amplitude + .5;

%%

% figure(7)
im = noise + grat;
im(~mask) = nan; 
% imagesc(im);


%  min(min(im )) * 255
%  max(max(im)) * 255
michelson_con = (max(max(im)) - min(min(im)))/( min(min(im)) + max(max(im)));


% specn = fftshift(fft2(noise));
% 
% figure(3); clf;
% imagesc(freq, freq, pow2db(abs(specn).^2));
% caxis(ca);
% title('power in noise');
% 
% 
% %% signal to noise ratio
% snr = abs(spec).^2 / abs(specn).^2;
% % 
% % figure(4); clf; 
% % imagesc(freq, freq, pow2db());
% 
% %% radial profiles
% pw = RadialProfile(abs(spec).^2, imsize/2, 30); %why imsize/2? Number of samples.
% pwn = RadialProfile(abs(specn).^2, imsize/2, 30);
% 
% snr = pw ./ pwn;
%     
% freq = 60 / pixelAngle * (0:imsize/2-1)/imsize;
% 
% snr_stim = pow2db(interp1(freq, snr, SpatialFreq, 'linear'));
% disp(snr_stim);
% 
% figure(5); clf; hold on;
% plot(freq, pow2db(pw), 'b-');
% plot(freq, pow2db(pwn), 'k-');
% xlabel('spatial frequency (cpd)'); ylabel('power (db)');
% set(gca, 'XScale', 'log', 'XTick', [1, 2, 10, 30]);
% 
% figure(6); clf; 
% plot(freq, pow2db(snr));
% xlabel('spatial frequency (cpd)'); ylabel('SNR (db)');
% set(gca, 'XScale', 'log', 'XTick', [1, 2, 10, 30]);
% end