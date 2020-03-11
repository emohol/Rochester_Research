path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\TIME_EXP_DATA';
% subject = 'A013';
% subject = 'Nikunj';
subject = 'A092';
% subject = 'A036';
pathtodata = fullfile(path,subject);

fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures_TIME';
sub_fig_path = fullfile(fig_path,subject,'psychfit');
% sub_fig_path = fullfile(fig_path,subject,'land30');
if (~isfolder(sub_fig_path))
    mkdir(sub_fig_path)
end


read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

[valid,counter] = countingTrialsNK(ppt);

%%
bin = 0; %set to 1 if binning
trials = struct;
count = 1;
for i=1:length(ppt)
    if valid.validTrials(i)
%     if valid.drift(i)
        trials(count).contrast = ppt{i,1}.contrast;
%         trials(count).contrast = Untitled3(ppt{i,1}.contrast,ppt{i,1}.eccentricity,ppt{i,1}.spatialFreq);
        trials(count).resp = ppt{i,1}.resp;
        trials(count).present = ppt{i,1}.present;
        trials(count).respTime = ppt{i,1}.responseTime;
        trials(count).ecc = ppt{i,1}.eccentricity;
        trials(count).presTime = ppt{i,1}.presTime;
        trials(count).spFreq = ppt{i,1}.spatialFreq;

        trials(count).pixelAngle = ppt{i,1}.pixelAngle;
        trials(count).xoffset = ppt{i,1}.xoffset;
        trials(count).yoffset = ppt{i,1}.yoffset;
        trials(count).cuePrf = ppt{i,1}.cuePrf;
        trials(count).fixOn = ppt{i,1}.fixOn;
        trials(count).presTime = ppt{i,1}.presTime;

        trials(count).cueCtr = ppt{i,1}.cueCtr;
        trials(count).saccOn = ppt{i,1}.saccOn;
        trials(count).saccOff = ppt{i,1}.saccOff;
        trials(count).flashOn = ppt{i,1}.flashOn;
        trials(count).rampOff = ppt{i,1}.rampOff;
        trials(count).stimOff = ppt{i,1}.stimOff;
        trials(count).quit = ppt{i,1}.quit; 
        
        %Correct for erroenous PEST levels
        if trials(count).contrast == 0 && trials(count).present == 1
            trials(count).present = 0;
        end
%         if trials(count).contrast ~= 0 && trials(count).present == 0
%             trials(count).contrast = 0;
%         end
        
        if trials(count).contrast > 0.5
            trials(count).contrast = 0.5;
        end
        count = count+1;
    end
end

%round off
if (bin==1)   
    for i=1:length(trials)
        trials(i).contrast = round([trials(i).contrast],2);
    end
else
    bin = 0;
end
ecc_levels = unique([trials.ecc]);
contrast = [trials.contrast];
contrast_levs = unique(contrast);
pres_Time = [trials.presTime];
sp_Freq = [trials.spFreq];
% falseAlarms = sum([trials.present]==0 &[trials.resp]==1);
% lures = sum([trials.present]==0);
% x=sort([trials.contrast]);

count = 1;
for i=1:length(ecc_levels)
    In1=([trials.ecc]==ecc_levels(i)); %indices%    
%     pres_levels= unique(pres_Time(In1));
    spFreq_lev = unique(sp_Freq(In1));
    [~,t] = size(spFreq_lev);
%     for p = 1:length(pres_levels) %different presentation times   
    for j = 1:t %different spatial frequencies
%         In2=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)); %indices
        In2=([trials.spFreq]==spFreq_lev(j) & [trials.ecc]==ecc_levels(i)); %indices
%         freq_levels = unique(sp_Freq(In2));
%         for c=1:length(freq_levels)
            plotting(count).ecc = ecc_levels(i);
%             plotting(count).presTime = pres_levels(p);
            plotting(count).spFreq = spFreq_lev(j);
%             In3=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1); %indices
            In3_present=([trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1); %indices
%             con_lev= unique(contrast(In3)); 
            pres_levels= unique(pres_Time(In3_present)); 
            plotting(count).presTime = pres_levels;
            x_present = find(In3_present); %Convert Logical into Sequential Index Vector
            plotting(count).contrast = trials(x_present(end)).contrast; %ALL TRIALS IN A CONDITION HAVE SAME CONTRAST            
%             for k = 1:length(con_lev)
            for k = 1:length(pres_levels)
                noTrials = sum([trials.presTime]==pres_levels(k) & [trials.ecc]==ecc_levels(i)& [trials.spFreq]==spFreq_lev(j));
                noStim = sum([trials.present]==1 & [trials.presTime]==pres_levels(k) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j));%stimulus displayed
                hits = sum([trials.present]==1 &[trials.resp]==1 & [trials.presTime]==pres_levels(k) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j));
                misses = sum([trials.present]==1 &[trials.resp]==0 & [trials.presTime]==pres_levels(k) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j));
                
                falseAlarms = sum([trials.present]==0 &[trials.resp]==1 & [trials.presTime]==pres_levels(k));
                lures = sum([trials.present]==0 & [trials.presTime]==pres_levels(k));
                
                plotting(count).noOfTrials(k) = noTrials; 
                plotting(count).noOfStimDisp(k) = noStim;
                plotting(count).hits(k) = hits;
                plotting(count).misses(k) = misses;
                plotting(count).percentCorrect(k) = 100 * (hits/noStim) ;
                
                plotting(count).falseAlarms(k) = falseAlarms;
                plotting(count).lures(k) = lures;                
                [plotting(count).dprime(k), ~] = dprime(hits,falseAlarms,noStim,lures);                
            end  
            count = count+1;
%         end
    end
end
% lvl = plotting(1).contrast; % contrast
%      hits = double(plotting(1).hits); % # correct
%      ulvl = unique(lvl); % unique contrast levels
%      nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
%      ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
%             
%     fh=1;
%     gi=1;
%      % do the psychometric fitting for these trials - get contrast
%      % thresholds
%      [thresh, par] = psyfit(lvl, hits,...
%          'Thresh', .75, 'Extra',...
%          'PlotHandle', fh+gi, 'Lapses', 'Auto');
%      title(sprintf('Ecc 0'));
%      xlim([0, 1]);
% 
%      x = logspace(-5, -1, 100); % range of contrast values to plot over
%      p = psyfun(x, par(1), par(2), 0.5, par(3), false, false); % use parameters from fit to plot curve
% 
%      figure(fh);
%      scatter(ulvl, ucorr, 10+5* nlvl, 'filled'); % plot data points sized by number of trials
%      hax(gi) = plot(x, p, '-', 'linewidth', 2); % plot the best fit psychometric curve
%      xlim([0, .1]);

%% FA vs Pres Time
fh=1;
figure(fh)
pres_levels_all = unique(pres_Time);
for k = 1:length(pres_levels_all)    
    falseAlarms = sum([trials.present]==0 &[trials.resp]==1 & [trials.presTime]==pres_levels_all(k));
    lures = sum([trials.present]==0 & [trials.presTime]==pres_levels_all(k));
    FA.falseAlarms(k) = falseAlarms;
    FA.lures(k) = lures;
end
% scatter(pres_levels_all,FA.falseAlarms./FA.lures,'filled')
xlabel('Pres. Time(ms)')
ylabel('FA Rate')
title('FA Rate vs Presentation Time (ALL Sessions)')
set(gca,'box','off','tickdir','out','FontSize',20)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%draw horizontal mean
l1=line([pres_levels_all(1) pres_levels_all(end)],[nanmean(FA.falseAlarms./FA.lures) nanmean(FA.falseAlarms./FA.lures)]);%,...
%     'DisplayName',['Mean - ',num2str(nanmean(FA.falseAlarms./FA.lures))]);
text(double(pres_levels_all(end)) , nanmean(FA.falseAlarms./FA.lures), ...
        ['Mean - ',num2str(nanmean(FA.falseAlarms./FA.lures))],'Color',[1 .0 .0],...
        'HorizontalAlignment','Center');
hold on
% legend(l1)
%draws circles for displaying data points below
for ii=1:length(pres_levels_all)
    plot( pres_levels_all(ii) , FA.falseAlarms(ii)/FA.lures(ii),...
        'o','MarkerSize',20,'Color',...
        [.3 .8 .7]);%,'LineWidth',ceil(FA.lures(ii)/10));
    %SETS NUMBERS INSIDE DATA POINTS NK
    text(double(pres_levels_all(ii)) , FA.falseAlarms(ii)/FA.lures(ii), ...
        num2str(FA.lures(ii)),'Color',[1 .0 .0],...
        'HorizontalAlignment','Center');
    hold on
end

%FA vs Pres Time per session
sess=[1];
count = 1;
if (read_Data == 1)
    for i=1:numel(files)
        if i==numel(files)
            break
        end
        if ~(strcmp(files(i).name(1:13),files(i+1).name(1:13)))
            sess = [sess; i+1];
            count = count+1;
        end
    end
else    
    for i=1:numel(filenames)
        if i==numel(filenames)
            break
        end
        if ~(strcmp(filenames(i).name(1:13),filenames(i+1).name(1:13)))
            sess = [sess; i+1];
            count = count+1;
        end
    end
end
saveas(figure(fh),[sub_fig_path '\FA_Pres_ALL'],'epsc');

hold off
fh=fh+1;
figure(fh)
for j=1:length(sess)
    subplot(length(sess),1,j)
    hold on
    if j==length(sess)
        pres_levels_all = unique(pres_Time(sess(j):length(trials)));
        for k = 1:length(pres_levels_all)
            falseAlarms = sum([trials(sess(j):length(trials)).present]==0 &[trials(sess(j):length(trials)).resp]==1 & [trials(sess(j):length(trials)).presTime]==pres_levels_all(k));
            lures = sum([trials(sess(j):length(trials)).present]==0 & [trials(sess(j):length(trials)).presTime]==pres_levels_all(k));
            FA_S(j).falseAlarms(k) = falseAlarms;
            FA_S(j).lures(k) = lures;
        end
    else
        pres_levels_all = unique(pres_Time(sess(j):sess(j+1)-1));
        for k = 1:length(pres_levels_all)
            falseAlarms = sum([trials(sess(j):sess(j+1)-1).present]==0 &[trials(sess(j):sess(j+1)-1).resp]==1 & [trials(sess(j):sess(j+1)-1).presTime]==pres_levels_all(k));
            lures = sum([trials(sess(j):sess(j+1)-1).present]==0 & [trials(sess(j):sess(j+1)-1).presTime]==pres_levels_all(k));
            FA_S(j).falseAlarms(k) = falseAlarms;
            FA_S(j).lures(k) = lures;
        end
    end
    
    % scatter(pres_levels_all,FA.falseAlarms./FA.lures,'filled')
    %draw horizontal mean
    line([pres_levels_all(1) pres_levels_all(end)],[nanmean(FA_S(j).falseAlarms./FA_S(j).lures) nanmean(FA_S(j).falseAlarms./FA_S(j).lures)])
    text(double(pres_levels_all(end)) , nanmean(FA_S(j).falseAlarms./FA_S(j).lures), ...
        ['Mean - ',num2str(nanmean(FA_S(j).falseAlarms./FA_S(j).lures))],'Color',[1 .0 .0],...
        'HorizontalAlignment','Center');
    xlabel('Pres. Time(ms)')
    ylabel('FA Rate')
    title(['FA Rate vs Presentation Time - Session ',num2str(j)])
    set(gca,'box','off','tickdir','out','FontSize',20)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

    hold on
    %draws circles for displaying data points below
    for ii=1:length(pres_levels_all)
        plot( pres_levels_all(ii) , FA_S(j).falseAlarms(ii)/FA_S(j).lures(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7]);%,'LineWidth',ceil(FA.lures(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(double(pres_levels_all(ii)) , FA_S(j).falseAlarms(ii)/FA_S(j).lures(ii), ...
            num2str(FA_S(j).lures(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
        hold on
    end
end
saveas(figure(fh),[sub_fig_path '\FA_Pres_Sess'],'epsc');


%%
hold off
fh=fh+1;    
% for i=1:length(plotting)
ind = [1:2:5];
for i=1:length(ind)
    figure(fh)
    subplot(3,2,2*i-1)
    hold on
    scatter(plotting(ind(i)).presTime,plotting(ind(i)).hits./plotting(ind(i)).noOfStimDisp,'filled')
    xlabel('Pres. Time(ms)')
    ylabel('Prop. Reported Present')
    set(gca,'box','off','tickdir','out','FontSize',20)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(['Ecc. - ',num2str(plotting(ind(i)).ecc),',Freq. - ',num2str(plotting(ind(i)).spFreq),',Con. - ', num2str(plotting(ind(i)).contrast)])
    
    subplot(3,2,2*i)
    hold on
    for ii=1:length(plotting(ind(i)).presTime)
        plot( double(plotting(ind(i)).presTime(ii)) , plotting(ind(i)).hits(ii)/plotting(ind(i)).noOfStimDisp(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7]);%,'LineWidth',ceil(FA.lures(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(double(plotting(ind(i)).presTime(ii)) ,  plotting(ind(i)).hits(ii)/plotting(ind(i)).noOfStimDisp(ii), ...
            num2str(plotting(ind(i)).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
    end
%     hold on
end
saveas(figure(fh),[sub_fig_path '\Per_Yes_SpFreq_2'],'epsc');

hold off
fh=fh+1;    
% for i=1:length(plotting)
ind = [2:2:6];
for i=1:length(ind)
    figure(fh)
    subplot(3,2,2*i-1)
    hold on
    scatter(plotting(ind(i)).presTime,plotting(ind(i)).hits./plotting(ind(i)).noOfStimDisp,'filled')
    xlabel('Pres. Time(ms)')
    ylabel('Prop. Reported Present')
    set(gca,'box','off','tickdir','out','FontSize',20)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(['Ecc. - ',num2str(plotting(ind(i)).ecc),',Freq. - ',num2str(plotting(ind(i)).spFreq),',Con. - ', num2str(plotting(ind(i)).contrast)])
    
    
    subplot(3,2,2*i)
    hold on
    for ii=1:length(plotting(ind(i)).presTime)
        plot( double(plotting(ind(i)).presTime(ii)) , plotting(ind(i)).hits(ii)/plotting(ind(i)).noOfStimDisp(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7]);%,'LineWidth',ceil(FA.lures(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(double(plotting(ind(i)).presTime(ii)) ,  plotting(ind(i)).hits(ii)/plotting(ind(i)).noOfStimDisp(ii), ...
            num2str(plotting(ind(i)).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
    end
end
saveas(figure(fh),[sub_fig_path '\Per_Yes_SpFreq_10'],'epsc');
%%
%draw Spatial Freq-wise, put the threshold value in the legend for now. Get
%the handle from the psyfit to adjust the color, legend, etc.
hold off
fh = fh+1;
ind = [1:2:5];
for i=1:length(ind)
    figure(fh)
    lvl = plotting(ind(i)).presTime; % contrast
    hits = plotting(ind(i)).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate=nanmean(plotting(ind(i)).falseAlarms./plotting(ind(i)).lures);
    subplot(3,1,i)
    hold on
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    [thresh, par]=psyfit_time(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ' ,Freq. - ',num2str(plotting(ind(i)).spFreq), ',Con. - ', num2str(plotting(ind(i)).contrast)],'PlotHandle',fh,'Extra',...
        'Chance',fa_rate,'Lapses','Auto','Thresh',0.75);
    plotting(ind(i)).thresh = round(thresh,2);%Round because that's the resolution of the shader
     xlabel('Pres. Time(ms)')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\SpFreq2'],'epsc');

hold off
fh = fh+1;
ind = [2:2:6];
for i=1:length(ind)
    figure(fh)
    lvl = plotting(ind(i)).presTime; % contrast
    hits = plotting(ind(i)).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate=nanmean(plotting(ind(i)).falseAlarms./plotting(ind(i)).lures);
    subplot(3,1,i)
    hold on
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    [thresh, par]=psyfit_time(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ' ,Freq. - ',num2str(plotting(ind(i)).spFreq), ',Con. - ', num2str(plotting(ind(i)).contrast)],'PlotHandle',fh,'Extra',...
        'Chance',fa_rate,'Lapses','Auto','Thresh',0.75);
    plotting(ind(i)).thresh = round(thresh,2);%Round because that's the resolution of the shader
     xlabel('Pres. Time(ms)')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\SpFreq10'],'epsc');

% hold off
% fh = fh+1;
% for i=5:6
%     figure(fh)
%     lvl = plotting(i).presTime; % contrast
%     hits = plotting(i).hits; % # correct
% %     tr = plotting(i).noOfTrials;
%     tr = plotting(i).noOfStimDisp;
%     fa_rate=max(plotting(i).falseAlarms./plotting(i).lures);
%     subplot(2,1,(i-4))
%     hold on
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Freq. - ',num2str(plotting(i).spFreq), ',Con. - ', num2str(plotting(i).contrast)],'PlotHandle',fh,'Extra',...
%         'Chance',fa_rate,'Lapses','Auto','Thresh',0.75);
%     plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader
% %      xlim([0, 0.5]);
% %      xlabel('Contrast')
%      xlabel('Pres. Time(ms)')
%      ylabel('Proportion Reported Present')
% %      title(subject)
% %      hold on
% %      scatter(lvl, hits, 10+5* tr);
% %     hold off
%     ulvl = unique(lvl); % unique contrast levels
%     nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
% %     nlvl = plotting(i).noOfTrials;
%     ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
%     set(gca,'box','off','tickdir','out','FontSize',20)
% %     hold on
% %     if (i==8)
% %         hold off
% %     end
% end
% % set(gca,'box','off','tickdir','out')
% % Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% saveas(figure(fh),[sub_fig_path '\Ecc8'],'epsc');


%% D-prime
fh = 4;
for i=1:2
    figure(fh)
    lvl = plotting(i).presTime; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    subplot(2,1,i)
    plot(plotting(i).presTime, plotting(i).dprime,'LineWidth',2,'Color',[0 0.5 0.5])
    hold on
    %draws circles for displaying data points below
    for ii=1:length(plotting(i).presTime)
        plot( plotting(i).presTime(ii) , plotting(i).dprime(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7],'LineWidth',ceil(plotting(i).noOfStimDisp(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(double(plotting(i).presTime(ii)) , plotting(i).dprime(ii), ...
            num2str(plotting(i).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
        hold on
    end
    %DISPLAY THRESHOLD AND TOTAL NUMBER OF TRIALS WITH STIMULUS DISPLAYED
    fontsize = 14;
    text( double(plotting(i).presTime(end)) , 0,...
        strcat('threshold =',num2str(round(plotting(i).thresh*100)/100)),...
        'FontWeight','Bold','FontSize',fontsize,'Color',[.3 .3 .5]);
    text( double(plotting(i).presTime(end)) , -0.5,...
        strcat('N = ',num2str(sum(plotting(i).noOfStimDisp))),...
        'FontSize',fontsize,'Color',[.3 .3 .5]);
    
    % Here we draw a line to mark the threshold
    [~,y2] = min(abs(plotting(i).presTime - plotting(i).thresh));
    %Vertical line
    line( [plotting(i).thresh plotting(i).thresh] , [min([plotting(i).dprime 0]) plotting(i).dprime(y2)], ...
        'color','r','LineStyle',':', 'LineWidth',2);
    %Horizontal line
    line( [0 plotting(i).thresh] , [plotting(i).dprime(y2) plotting(i).dprime(y2)],...
        'color','r','LineStyle',':', 'LineWidth',2);
    
%      xlim([0, 0.5]);
     xlabel('Pres. Time(ms)')
     ylabel('D-Prime')
     title(['Ecc. - ',num2str(plotting(i).ecc)',',Freq. - ',num2str(plotting(i).spFreq),',Con. - ',num2str(plotting(i).contrast),])
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc0_Dprime'],'epsc');


fh = fh+1;
for i=3:4
    figure(fh)
    lvl = plotting(i).presTime; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    subplot(2,1,i-2)
    plot(plotting(i).presTime, plotting(i).dprime,'LineWidth',2,'Color',[0 0.5 0.5])
    hold on
    %draws circles for displaying data points below
    for ii=1:length(plotting(i).presTime)
        plot( plotting(i).presTime(ii) , plotting(i).dprime(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7],'LineWidth',ceil(plotting(i).noOfStimDisp(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(double(plotting(i).presTime(ii)) , plotting(i).dprime(ii), ...
            num2str(plotting(i).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
        hold on
    end
    %DISPLAY THRESHOLD AND TOTAL NUMBER OF TRIALS WITH STIMULUS DISPLAYED
    fontsize = 14;
    text( double(plotting(i).presTime(end)) , 0,...
        strcat('threshold =',num2str(round(plotting(i).thresh*100)/100)),...
        'FontWeight','Bold','FontSize',fontsize,'Color',[.3 .3 .5]);
    text( double(plotting(i).presTime(end)) , -0.5,...
        strcat('N = ',num2str(sum(plotting(i).noOfStimDisp))),...
        'FontSize',fontsize,'Color',[.3 .3 .5]);
    
    % Here we draw a line to mark the threshold
    [~,y2] = min(abs(plotting(i).presTime - plotting(i).thresh));
    %Vertical line
    line( [plotting(i).thresh plotting(i).thresh] , [min([plotting(i).dprime 0]) plotting(i).dprime(y2)], ...
        'color','r','LineStyle',':', 'LineWidth',2);
    %Horizontal line
    line( [0 plotting(i).thresh] , [plotting(i).dprime(y2) plotting(i).dprime(y2)],...
        'color','r','LineStyle',':', 'LineWidth',2);
    
%      xlim([0, 0.5]);
     xlabel('Pres. Time (ms)')
     ylabel('D-Prime')
     title(['Ecc. - ',num2str(plotting(i).ecc)',',Freq. - ',num2str(plotting(i).spFreq),',Con. - ',num2str(plotting(i).contrast),])
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc4_Dprime'],'epsc');

fh = fh+1;
for i=5:6
    figure(fh)
    lvl = plotting(i).presTime; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    subplot(2,1,i-4)
    plot(plotting(i).presTime, plotting(i).dprime,'LineWidth',2,'Color',[0 0.5 0.5])
    hold on
    %draws circles for displaying data points below
    for ii=1:length(plotting(i).presTime)
        plot( plotting(i).presTime(ii) , plotting(i).dprime(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7],'LineWidth',ceil(plotting(i).noOfStimDisp(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(double(plotting(i).presTime(ii)) , plotting(i).dprime(ii), ...
            num2str(plotting(i).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
        hold on
    end
    %DISPLAY THRESHOLD AND TOTAL NUMBER OF TRIALS WITH STIMULUS DISPLAYED
    fontsize = 14;
    text(double(plotting(i).presTime(end)) , 0,...
        strcat('threshold =',num2str(round(plotting(i).thresh*100)/100)),...
        'FontWeight','Bold','FontSize',fontsize,'Color',[.3 .3 .5]);
    text( double(plotting(i).presTime(end)) , -0.5,...
        strcat('N = ',num2str(sum(plotting(i).noOfStimDisp))),...
        'FontSize',fontsize,'Color',[.3 .3 .5]);
    
    % Here we draw a line to mark the threshold
    [~,y2] = min(abs(plotting(i).presTime - plotting(i).thresh));
    %Vertical line
    line( [plotting(i).thresh plotting(i).thresh] , [min([plotting(i).dprime 0]) plotting(i).dprime(y2)], ...
        'color','r','LineStyle',':', 'LineWidth',2);
    %Horizontal line
    line( [0 plotting(i).thresh] , [plotting(i).dprime(y2) plotting(i).dprime(y2)],...
        'color','r','LineStyle',':', 'LineWidth',2);
    
%      xlim([0, 0.5]);
     xlabel('Pres. Time (ms)')
     ylabel('D-Prime')
     title(['Ecc. - ',num2str(plotting(i).ecc)',',Freq. - ',num2str(plotting(i).spFreq),',Con. - ',num2str(plotting(i).contrast),])
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc8_Dprime'],'epsc');

%% Sensitivity/Threshold summary plots
%Pres time vs sensitivity (2 cpd)
fh=fh+1;
figure(fh)
% for c=1:length(freq_levels)
%     for p = 1:length(pres_levels)
        for i=1:length(ecc_levels)
            In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
            In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
            t1 = plotting(In1).thresh;
            t2 = plotting(In2).thresh;
            if t1>0.5
                t1 = 0.5;
            elseif t1<0
                t1 = 0;
            end
            if t2>0.5
                t2 = 0.5;
            elseif t2<0
                t2 = 0;
            end
            sens = [1/t1 1/t2];
%             if sens>0.5
%                 sens = 0.5;
%             elseif sens<0
%                 sens = 0;
%             end
            if ecc_levels(i)==0
                plot([50 500],sens,'-s','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
                hold on
            elseif ecc_levels(i)==4
                plot([50 500],sens,'-v','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
                hold on
            elseif ecc_levels(i)==8
                plot([50 500],sens,'-o','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
                hold on
%             elseif ecc_levels(i)==0 && freq_levels(c)==10
%                 plot([50 500],th,'--s','linewidth', 2,'Color','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==4 && freq_levels(c)==10
%                 plot([50 500],th,'--v','linewidth', 2,'Color','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==8 && freq_levels(c)==10
%                 plot([50 500],th,'--o','linewidth', 2,'Color','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
            end
        end
%     end
    
% end
hold off
set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend
xticks([50 500])
xlim([-50 600])
xlabel('Presentation time')
ylabel('Sensitivity')
title('Low SpFreq')
saveas(figure(fh),[sub_fig_path '\Sens_Pres_2'],'epsc');
%%
%Pres time vs sensitivity (10 cpd)

alpha = [];
fh=fh+1;
figure(fh)
% y_tick=[];
% for c=1:length(freq_levels)
%     for p = 1:length(pres_levels)
        for i=1:length(ecc_levels)
            In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
            In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
            t1 = plotting(In1).thresh;
            t2 = plotting(In2).thresh;
            if t1>0.5
                t1 = 0.5;
            elseif t1<0
                t1 = 0;
            end
            if t2>0.5
                t2 = 0.5;
            elseif t2<0
                t2 = 0;
            end
            sens = [1/t1 1/t2];
%             alpha = [alpha atand((sens(2) - sens(1))/(450))];            
            alpha = [alpha (sens(2) - sens(1))/(450)];
            formatSpec = "(sens/ms)";
%             if ecc_levels(i)==0 && freq_levels(c)==2
%                 plot([50 500],th,'-s','linewidth', 2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==4 && freq_levels(c)==2
%                 plot([50 500],th,'-v','linewidth', 2,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==8 && freq_levels(c)==2
%                 plot([50 500],th,'-o','linewidth', 2,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
            if ecc_levels(i)==0
                plot([50 500],sens,'-s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 txt = sprintf(formatSpec,alpha(i));
                txt = ['\angle' alpha(i) formatSpec];
%                 txt = ['\angle',alpha(i)];
                text(50,sens(1),txt,'FontSize',20)
                hold on
            elseif ecc_levels(i)==4
                plot([50 500],sens,'-v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
                txt = ['\angle' alpha(i) formatSpec];
                text(50,sens(1),txt,'FontSize',20)
                hold on
            elseif ecc_levels(i)==8
                plot([50 500],sens,'-o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
                txt = ['\angle' alpha(i) formatSpec];
                text(50,sens(1),txt,'FontSize',20)
                hold on
            end
%             y_tick = [y_tick th];
        end
%     end
%     
% end
hold off
set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend
xticks([50 500])
xlim([-50 600])
% ylim([0 0.5])
% yticks(sort(y_tick))
xlabel('Presentation time')
ylabel('Sensitivity')
title('High SpFreq')
saveas(figure(fh),[sub_fig_path '\Sens_Pres_10'],'epsc');
%%

%Pres time vs threshold (2 cpd)
fh=fh+1;
figure(fh)
% for c=1:length(freq_levels)
%     for p = 1:length(pres_levels)
        for i=1:length(ecc_levels)
            In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
            In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
            t1 = plotting(In1).thresh;
            t2 = plotting(In2).thresh;
            if t1>0.5
                t1 = 0.5;
            elseif t1<0
                t1 = 0;
            end
            if t2>0.5
                t2 = 0.5;
            elseif t2<0
                t2 = 0;
            end
            thresh = [t1 t2];
            if ecc_levels(i)==0
                plot([50 500],thresh,'-s','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
                hold on
            elseif ecc_levels(i)==4
                plot([50 500],thresh,'-v','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
                hold on
            elseif ecc_levels(i)==8
                plot([50 500],thresh,'-o','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
                hold on
%             elseif ecc_levels(i)==0 && freq_levels(c)==10
%                 plot([50 500],th,'--s','linewidth', 2,'Color','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==4 && freq_levels(c)==10
%                 plot([50 500],th,'--v','linewidth', 2,'Color','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==8 && freq_levels(c)==10
%                 plot([50 500],th,'--o','linewidth', 2,'Color','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
            end
        end
%     end
    
% end
hold off
set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend
xticks([50 500])
xlim([-50 600])
ylim([0 0.5])
xlabel('Presentation time')
ylabel('Threshold')
title('Low SpFreq')
saveas(figure(fh),[sub_fig_path '\Thresh_Pres_2'],'epsc');

%Pres time vs threshold (10 cpd)
fh=fh+1;
figure(fh)
% y_tick=[];
% for c=1:length(freq_levels)
%     for p = 1:length(pres_levels)
        for i=1:length(ecc_levels)
            In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
            In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
            t1 = plotting(In1).thresh;
            t2 = plotting(In2).thresh;
            if t1>0.5
                t1 = 0.5;
            elseif t1<0
                t1 = 0;
            end
            if t2>0.5
                t2 = 0.5;
            elseif t2<0
                t2 = 0;
            end
            thresh = [t1 t2];
            
%             if ecc_levels(i)==0 && freq_levels(c)==2
%                 plot([50 500],th,'-s','linewidth', 2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==4 && freq_levels(c)==2
%                 plot([50 500],th,'-v','linewidth', 2,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
%             elseif ecc_levels(i)==8 && freq_levels(c)==2
%                 plot([50 500],th,'-o','linewidth', 2,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
%                 hold on
            if ecc_levels(i)==0
                plot([50 500],thresh,'-s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
                hold on
            elseif ecc_levels(i)==4
                plot([50 500],thresh,'-v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
                hold on
            elseif ecc_levels(i)==8
                plot([50 500],thresh,'-o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
                    'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
                hold on
            end
%             y_tick = [y_tick th];
        end
%     end
%     
% end
hold off
set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend
xticks([50 500])
xlim([-50 600])
ylim([0 0.5])
% yticks(sort(y_tick))
xlabel('Presentation time')
ylabel('Threshold')
title('High SpFreq')
saveas(figure(fh),[sub_fig_path '\Thresh_Pres_10'],'epsc');
%% Radial Colorbar plot of thresholds

imsize_x = 1080;
imsize_y = 1920;
pixelAngle = 1.06; % arcmin per pixel, was 0.8 initially
width = 1;
fh=fh+1;
fig_count =fh;
for p = 1:length(pres_levels) %different presentation times
    for c=1:length(freq_levels)
        grat = zeros(imsize_x,imsize_y);        
        for ii = 1:imsize_x
            for jj = 1:imsize_y
                r = sqrt ((ii - imsize_x/2)^2 + (jj - imsize_y/2)^2)* pixelAngle / 60;
                if ((r >= ecc_levels(1) - width/2) && (r <= ecc_levels(1) + width/2))
                    In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(1) & [plotting.spFreq]==freq_levels(c)); %indices
                    scale = plotting(In1).thresh;
                elseif ((r >= ecc_levels(2) - width/2) && (r <= ecc_levels(2) + width/2))
                    In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(2) & [plotting.spFreq]==freq_levels(c)); %indices
                    scale = plotting(In1).thresh;
                elseif ((r >= ecc_levels(3) - width/2) && (r <= ecc_levels(3) + width/2))
                    In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(3) & [plotting.spFreq]==freq_levels(c)); %indices
                    scale = plotting(In1).thresh;
                else
                    scale = 0;
                end
                grat(ii, jj) = scale;
            end
        end
        figure(fig_count)
        imagesc(grat);
        set (gca, 'CLim',[0 0.5])
        colorbar
        title(['Pres. - ',num2str(pres_levels(p)),' ,Freq. - ',num2str(freq_levels(c))])
        saveas(figure(fig_count),[sub_fig_path '\Color_' ['Pres_',num2str(pres_levels(p)),'_Freq_',num2str(freq_levels(c))]],'epsc');
        fig_count =fig_count+1;
    end
end




%%
%Save threshold summary table

% thresh_T = table;
% formatSpec = 'Ecc-%d, Pres-%d, SpFreq-%d';
% dummy_cell = {};
% for i =1:12
%     name = sprintf(formatSpec,plotting(i).ecc,plotting(i).presTime,plotting(i).spFreq);
%     row_1 = {name,num2str(round(plotting(i).thresh,2))};
%     dummy_cell{i,1} = row_1;
%     colnames = {'Condition';'Threshold'};
%     dummy = cell2table(dummy_cell{i,1},'VariableNames',colnames);
%     thresh_T = [thresh_T; dummy];    
% end
% writetable(thresh_T,[sub_fig_path '\psychfit.csv'],'WriteVariableNames',true);

