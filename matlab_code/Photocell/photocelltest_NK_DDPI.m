%%
path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\BENQ\PhotoCellBenQ\with_eye_tracking\BlurRedONBlurBusterON_FINAL_DPI';
% path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\PhotoCellBenQ\BlurRedONBlurBusterON_FINAL_DDPI_LABJACK';
% subject = 'A013';
% subject = 'BlurRedPH5';
% pathtodata = fullfile(path,subject);

% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj_NoMask';
read_Data = 1; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(path, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(path, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

%%
% load('NoBlurBustersDPI.mat')
load('withBlurBusters.mat')
% load('noBlurBust.mat')
%%

Fs = 1024;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
% L = length(v);             % Length of signal
L = length(v(10000:end));             % Length of signal
t = (0:L-1)*T;        % Time vector
% t = (0:L-1)*T;        % Time vector



% Y = fft(v);
%DDPI
Y = fft(v(10000:end));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(1)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
set(gca,'box','off','tickdir','out','FontSize',30)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
title('Frequency spectrum - dDPI data')
% saveas(figure(1),'DPI_withBlurBusters_freq_3.jpg');
% saveas(figure(1),'dDPI_BlurBustersOff_freq_3.jpg');
% saveas(figure(1),'dDPI_withBlurBusters_freq_1.jpg');
% lum = data.photocell{end-100};

%DPI
% lum = v(1:end);
% t = (1:length(v)) / 1000;
% [~, locs] = findpeaks(double(v), 'minpeakheight', 1);

%FOR DDPI
% v=lum;
t = (10000:length(v)) / 1000;%Excluding first 10 secs
[~, locs] = findpeaks(double(v(10000:end)), 'minpeakheight', 1);




d = diff(t(locs));
mnT = mean(d);
fprintf('cycle frequency = %1.3fHz, refresh rate = %1.3fHz\n', ...
    1/mnT, 2/mnT);
fprintf('%i peaks in %1.3fs --> cycle frequency = %1.3fHz\n',...
    length(locs), t(end), length(locs) / t(end));

% miss_fr_high = sum(d>mnT);%Number of missed frames
% miss_fr_low = sum(d<mnT);%Number of missed frames

%200Hz range
% miss_fr_high = sum(d>0.0105);%Number of missed frames below 190Hz
% miss_fr_low = sum(d<0.0095);%Number of missed frames above 210Hz

%240Hz range
miss_fr_high = sum(d>0.0087);%Number of missed frames below 230Hz
miss_fr_low = sum(d<0.008);%Number of missed frames above 250Hz
miss_fr = miss_fr_high + miss_fr_low;


figure(2)
% plot(v)
plot(v(10000:end))
xlabel('Time (ms)');
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage');
set(gca,'box','off','tickdir','out','FontSize',30)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
title('dDPI data')
% saveas(figure(2),'DPI_raw_data_1.jpg');
% saveas(figure(2),'DPI_withBlurBusters_time_diff_3_raw_data.jpg');
% saveas(figure(2),'dDPI_BlurBustersOff_time_diff_3_raw_data.jpg');
% saveas(figure(2),'dDPI_withBlurBusters_time_diff_1_raw_data.jpg');

figure(3)
histogram(d)
% h1=histogram(d,'DisplayName',['Total missed frames(%)-',num2str(miss_fr*100/length(d)),'Time diff < mean-',num2str(miss_fr_low*100/length(d))]);
xlabel('Time between peaks');
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
ylabel('No. of Occurences');
hold on
vertLineThrough(mnT, 'k');
% vertLineThrough(0.0105, 'k');
% vertLineThrough(0.0095, 'k');
% line([1,length(d)],[mnT,mnT],'Color','red', 'LineWidth', 2)
formatSpec = "Mean - %1.5f";
txt = sprintf(formatSpec,mnT);
text(mnT,2500,txt,'FontSize',20);

% message = sprintf('Line #1\nThe second line.\nAnd finally a third line.');
% text(5, 5, message, 'FontSize', 20, 'Color', [.6 .2 .6]);

formatSpec = "Percent Total missed frames - %1.5f, \n Time difference less than mean - %1.5f,\n Time difference more than mean - %1.5f";
txt = sprintf(formatSpec,miss_fr*100/length(d),miss_fr_low*100/length(d), miss_fr_high*100/length(d));
text(max(d)/2,2500,txt,'FontSize',20);

% text(max(d), mnT, 'Mean - %d',mnT, 'FontSize', 20);
set(gca,'box','off','tickdir','out','FontSize',30)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% legend(h1)
title('Time difference between peaks')
% saveas(figure(3),'dDPI_withBlurBusters_time_diff_3.jpg');
% saveas(figure(3),'DPI_withBlurBusters_time_diff_3.jpg');
% saveas(figure(3),'dDPI_BlurBustersOff_time_diff_3.jpg');
% saveas(figure(3),'dDPI_withBlurBusters_time_diff_1.jpg');



