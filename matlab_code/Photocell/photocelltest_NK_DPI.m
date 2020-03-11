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
load('NoBlurBustersDPI.mat')
% load('withBlurBusters.mat')
% load('noBlurBust.mat')
% load('ddpi_asus_frameSkip.mat')
%%

Fs = 1024;            % Sampling frequency                    
T = 1/Fs;             % Sampling period    
L = length(v);             % Length of signal
t = (0:L-1)*T;        % Time vector
% t = (0:L-1)*T;        % Time vector



Y = fft(v);
P2 = abs(Y/L);
% Power = (abs(Y).^2)/L;
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
title('Amplitude Frequency spectrum')
saveas(figure(1),'dDPI_Asus_amp_3.jpg');
% saveas(figure(1),'DPI_BlurBustersOff_freq_3.jpg');
% lum = data.photocell{end-100};

figure(4)
plot(f,pow2db(P1.*P1)) 
title('Power Spectrum')
xlabel('f (Hz)')
ylabel('Power')
set(gca,'box','off','tickdir','out','FontSize',30)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
title('Power Frequency spectrum')
% saveas(figure(4),'dDPI_Asus_power_3.jpg');

%DPI
% lum = v(1:end);
t = (1:length(v)) / 1000;
[~, locs] = findpeaks(double(v), 'minpeakheight', 1);


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
miss_fr_high = sum(d>0.0111);%Number of missed frames below 180Hz
% miss_fr_low = sum(d<0.0095);%Number of missed frames above 210Hz
miss_fr_low = sum(d<0.0091);%Number of missed frames above 220Hz

%240Hz range
% miss_fr_high = sum(d>0.0087);%Number of missed frames below 230Hz
% miss_fr_low = sum(d<0.008);%Number of missed frames above 250Hz
miss_fr = miss_fr_high + miss_fr_low;


figure(2)
plot(v)
xlabel('Time (ms)');
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Voltage');
set(gca,'box','off','tickdir','out','FontSize',30)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
title('Raw data')
% saveas(figure(2),'dDPI_Asus_raw_data_3.jpg');
% saveas(figure(2),'DPI_withBlurBusters_time_diff_3_raw_data.jpg');
% saveas(figure(2),'DPI_BlurBustersOff_time_diff_3_raw_data.jpg');


figure(3)
histogram(d)
% h1=histogram(d,'DisplayName',['Total missed frames(%)-',num2str(miss_fr*100/length(d)),'Time diff < mean-',num2str(miss_fr_low*100/length(d))]);
xlabel('Time between peaks');
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 14, 'FontWeight', 'bold');
xlim([0 0.05])
ylabel('No. of Occurences');
hold on
vertLineThrough(mnT, 'k');
% vertLineThrough(0.0105, 'k');
% vertLineThrough(0.0095, 'k');
% line([1,length(d)],[mnT,mnT],'Color','red', 'LineWidth', 2)
hold on
formatSpec = "Mean - %1.5f";
txt = sprintf(formatSpec,mnT);
text(mnT,1500,txt,'FontSize',20);

% message = sprintf('Line #1\nThe second line.\nAnd finally a third line.');
% text(5, 5, message, 'FontSize', 20, 'Color', [.6 .2 .6]);
hold on
formatSpec = "Percent Total missed frames - %1.5f, \n Time difference less than mean - %1.5f,\n Time difference more than mean - %1.5f";
txt = sprintf(formatSpec,miss_fr*100/length(d),miss_fr_low*100/length(d), miss_fr_high*100/length(d));
% text(max(d)/2,1500,txt,'FontSize',20);
text(0.03,1500,txt,'FontSize',20);

% text(max(d), mnT, 'Mean - %d',mnT, 'FontSize', 20);
set(gca,'box','off','tickdir','out','FontSize',30)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% legend(h1)
title('Time difference between peaks')
% saveas(figure(3),'dDPI_Asus_time_diff_3.jpg');
% saveas(figure(3),'dDPI_Asus_time_diff_3.fig');
% saveas(figure(3),'DPI_withBlurBusters_time_diff_3.jpg');
% saveas(figure(3),'DPI_BlurBustersOff_time_diff_3.jpg');




