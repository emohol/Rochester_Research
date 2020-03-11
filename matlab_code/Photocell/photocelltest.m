%% photocell test
nfft = 512;
Fs = 1000;
allF = nan(length(data.photocell), nfft/2 + 1);
for ii = 1:length(data.photocell)
%     if counts.notracks(ii)
%         continue;
%     end
%     e = round(trials(ii).saccOn);
%     s = round(trials(ii).stimOff);
%     p = data.photocell{ii}(e:s);
    p = data.photocell{ii};
    
    [f, freq] = pwelch(p - mean(p), nfft, 128, nfft, 1000);
%     [f, freq] = pwelch(p - mean(p));
    allF(ii, :) = f;
    %     disp(freq(end));
    
    %     figure(1); clf;
    %     subplot(2, 1, 1);
    %     plot(p); xlim([500, 1000]);
    %     subplot(2, 1, 2);
    %     plot(freq, f);
    %     pause();
end

figure(10); clf;
boundedline(freq, nanmean(allF), nanstd(allF));
% [~, locs] = findpeaks(nanmean(allF), 'minpeakdistance', 10, 'Threshold', .5 * 10^-10);
% [~, locs] = findpeaks(nanmean(allF));
% hl = vertLineThrough(freq(locs));
% set(hl, 'linewidth', 1, 'color', 'k');
% xlim([0, 250]);
xlabel('frequency (hz)');
ylabel('amplitude');

%%
plot(lowpass(trials(1).photodiode,150,1000))

