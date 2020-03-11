% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\PhotoCellBenQ\BlurRedPH5';
path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\PhotoCellBenQ';
% subject = 'A013';
subject = 'Photocell3withBlurRed';
pathtodata = fullfile(path,subject);

% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj_NoMask';
read_Data = 1; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

% for i=1:length(ppt)
for i=1
% for i=1148:1198
    lum = data.photocell{1:end};
%     lum = data.photocell{end-100};
    t = (1:length(lum)) / 1000;

    [~, locs] = findpeaks(double(lum), 'minpeakheight', 1);



    d = diff(t(locs));
    mnT = mean(d);
    fprintf('cycle frequency = %1.3fHz, refresh rate = %1.3fHz\n', ...
        1/mnT, 2/mnT);
    fprintf('%i peaks in %1.3fs --> cycle frequency = %1.3fHz\n',...
        length(locs), t(end), length(locs) / t(end));

    figure(1); clf; hold on;
    plot(t*1000, lum, 'bo-', 'linewidth', 2);
    xlim(1000*t([1, end]));
    xlabel('time (ms)');
    ylabel('voltage');
    set(gcf, 'Color', 'w');
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');

    figure(2); clf; hold on;
    plot(t*1000, lum, 'b-', 'linewidth', 2);
    xlim(1000*t([1, end]));
    xlabel('time (ms)');
    set(gcf, 'Color', 'w');
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('voltage');
    vertLineThrough([13.67, 20.51], 'k');
    text(13.7, 1.8, '13.67ms', 'FontSize', 20);
    text(20.6, 1.8, '20.51ms', 'FontSize', 20);
%     name = sprintf('C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/test_new_photodiode/pd_%i.png',i);
    name = [pathtodata '\1.fig'];
    saveas(figure(2),name);    
end