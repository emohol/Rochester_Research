path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';

% subject = 'Nikunj';
% subject = 'A013';
% subject = 'A092';
subject = 'A036';
pathtodata = fullfile(path,subject);

fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures';
sub_fig_path = fullfile(fig_path,subject,'psychfit/perf50');
% sub_fig_path = fullfile(fig_path,subject,'psychfit');
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

nBoots = 1;
th_per = 0.5;%Threshold calculated at this performance

%%
bin = 1; %set to 1 if binning
trials = struct;
count = 1;
for i=1:length(ppt)
%     if valid.validTrials(i)
    if valid.drift(i)
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
        
        if trials(count).contrast > 0.6
            trials(count).contrast = 0.6;
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
pres_Time = [trials.presTime];
sp_Freq = [trials.spFreq];
% falseAlarms = sum([trials.present]==0 &[trials.resp]==1);
% lures = sum([trials.present]==0);
% x=sort([trials.contrast]);

count = 1;
for i=1:length(ecc_levels)
    In1=([trials.ecc]==ecc_levels(i)); %indices%    
    pres_levels= unique(pres_Time(In1));
    for p = 1:length(pres_levels) %different presentation times   
        In2=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)); %indices
        freq_levels = unique(sp_Freq(In2));
        for c=1:length(freq_levels)
            plotting(count).ecc = ecc_levels(i);
            plotting(count).presTime = pres_levels(p);
            plotting(count).spFreq = freq_levels(c);
            In3=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1); %indices
            con_lev= unique(contrast(In3)); 
            plotting(count).contrast = con_lev;            
            for j = 1:length(con_lev)
                noTrials = sum([trials.contrast]==plotting(count).contrast(j) & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)& [trials.spFreq]==freq_levels(c));
                noStim = sum([trials.contrast] == plotting(count).contrast(j) & [trials.present]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));%stimulus displayed
                hits = sum([trials.contrast]==plotting(count).contrast(j) & [trials.present]==1 &[trials.resp]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));
                misses = sum([trials.contrast]==plotting(count).contrast(j) & [trials.present]==1 &[trials.resp]==0 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));
                
                falseAlarms = sum([trials.present]==0 &[trials.resp]==1 & [trials.presTime]==pres_levels(p));
                lures = sum([trials.present]==0 & [trials.presTime]==pres_levels(p));
                
                plotting(count).noOfTrials(j) = noTrials; 
                plotting(count).noOfStimDisp(j) = noStim;
                plotting(count).hits(j) = hits;
                plotting(count).misses(j) = misses;
                plotting(count).percentCorrect(j) = 100 * (hits/noStim) ;
                
                plotting(count).falseAlarms(j) = falseAlarms;
                plotting(count).lures(j) = lures;                
                [plotting(count).dprime(j), ~] = dprime(hits,falseAlarms,noStim,lures);                
            end  
%             plotting(count).contrast(end+1)= 0;
%             plotting(count).noOfStimDisp(end+1)= sum([trials.present]==0 & [trials.presTime]==pres_levels(p));
%             plotting(count).hits(end+1)= sum([trials.present]==0 &[trials.resp]==0 & [trials.presTime]==pres_levels(p));
            count = count+1;
        end
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
%%
%draw eccentricity-wise, put the threshold value in the legend for now. Get
%the handle from the psyfit to adjust the color, legend, etc.


for i=1:4
    figure(1)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    fa_rate = plotting(i).falseAlarms(1)/plotting(i).lures(1);
    subplot(2,2,i)
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    [thresh, par, threshB, parB]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
        ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra',...
       'Chance', fa_rate,'Lapses','Auto','Thresh',th_per,'Boots', nBoots,'DistType','Normal');%,'Xlim',[0 0.6],'Ylim',[0 1]);
    plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader
%     plotting(i).thresh = round(mean(threshB),2);%Round because that's the resolution of the shader    
    plotting(i).std = round(std(threshB),2);%Round because that's the resolution of the shader
    plotting(i).par = par;
    plotting(i).threshB = threshB;
    plotting(i).parB = parB;
    ulvl = unique(lvl); % unique contrast levels
    nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
%     nlvl = plotting(i).noOfTrials;
    ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
%     hold on
%     if (i==4)
%         hold off
%     end
%     title(sprintf('Ecc 0'));
%      xlim([0, .1]);
     
%      x = logspace(-5, -1, 100); % range of contrast values to plot over
%      p = psyfun(x, par(1), par(2), 0.5, par(3), false, false); % use parameters from fit to plot curve
     
%      scatter(ulvl, ucorr, 10+5* nlvl); % plot data points sized by number of trials
%      hold on
%      scatter(lvl, hits, 10+5* tr); % plot data points sized by number of trials
%      plot(x, p, '-', 'linewidth', 2); % plot the best fit psychometric curve
     ylim([0, 1]);
%      xlabel('Contrast')
     xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
%      title(subject)
%      hold off
%     han(i).Color = 'blue';
%     han(i).DisplayName = ['Pres - ',num2str(plotting(i).presTime),',Ecc - ',num2str(plotting(i).ecc),...
%         ',SpFreq - ',num2str(plotting(i).spFreq),',Thresh - ',num2str(thresh)];
%     
%     pause
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

for i=5:8
    figure(2)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    fa_rate = plotting(i).falseAlarms(1)/plotting(i).lures(1);
    subplot(2,2,(i-4))
%     hold on
%     if i==6
%         [thresh, par, threshB, parB]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',2,'Extra',...
%         'Chance',fa_rate,'Lapses','Auto','Thresh',0.75,'Boots', nBoots,'DistType','Weibull','Xlim',[0 0.6],'Ylim',[0 1]);
%     else
        [thresh, par, threshB, parB]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
        ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',2,'Extra',...
        'Chance',fa_rate,'Lapses','Auto','Thresh',th_per,'Boots', nBoots,'DistType','Normal');%,'Xlim',[0 0.6],'Ylim',[0 1]);
%     end
    plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader
    plotting(i).std = round(std(threshB),2);%Round because that's the resolution of the shader
    plotting(i).par = par;
    plotting(i).threshB = threshB;
    plotting(i).parB = parB;
%      xlim([0, 0.6]);
%      xlabel('Contrast')
    ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported Present')
%      title(subject)
%     hold on
%     scatter(lvl, hits, 10+5* tr);
    ulvl = unique(lvl); % unique contrast levels
    nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
%     nlvl = plotting(i).noOfTrials;
    ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
    set(gca,'box','off','tickdir','out','FontSize',20)
%     hold on
%     if (i==8)
%         hold off
%     end
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

for i=9:12
    figure(3)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    fa_rate = plotting(i).falseAlarms(1)/plotting(i).lures(1);
    subplot(2,2,(i-8))
%     hold on
    [thresh, par, threshB, parB]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
        ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',3,'Extra',...
        'Chance',fa_rate,'Lapses','Auto','Thresh',th_per,'Boots', nBoots,'DistType','Normal');%,'Xlim',[0 0.6],'Ylim',[0 1]);
    plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader
    plotting(i).std = round(std(threshB),2);%Round because that's the resolution of the shader
    plotting(i).par = par;
    plotting(i).threshB = threshB;
    plotting(i).parB = parB;
%      xlim([0, 0.6]);
%      xlabel('Contrast')
    ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported Present')
%      title(subject)
%      hold on
%      scatter(lvl, hits, 10+5* tr);
%     hold off
    ulvl = unique(lvl); % unique contrast levels
    nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
%     nlvl = plotting(i).noOfTrials;
    ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
    set(gca,'box','off','tickdir','out','FontSize',20)
%     hold on
%     if (i==8)
%         hold off
%     end
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

%Make a summary table of thresholds

saveas(figure(1),[sub_fig_path '\Ecc0'],'epsc');
saveas(figure(2),[sub_fig_path '\Ecc4'],'epsc');
saveas(figure(3),[sub_fig_path '\Ecc8'],'epsc');
if nBoots>1
    matfilename = strcat(subject,'_plot_BOOT.mat');
    save(matfilename,'plotting')
else
    matfilename = strcat(subject,'_plot.mat');
    save(matfilename,'plotting')
end


%% Psychfit (Comparison across all eccentricities)
fh=4;
ind = [1:4:9];
for i=1:length(ind)
    figure(fh)
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
%     subplot(3,1,i)
    if i==length(ind)
        subplot(2,2,[3,4])
    else
        subplot(2,2,i)
    end
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');%,'Xlim',[0 0.6],'Ylim',[0 1]);
   ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Pres_50_spFreq_2'],'epsc');

fh=fh+1;
ind = [3:4:11];
for i=1:length(ind)
    figure(fh)
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
%     subplot(3,1,i)
    if i==length(ind)
        subplot(2,2,[3,4])
    else
        subplot(2,2,i)
    end
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');%,'Xlim',[0 0.6],'Ylim',[0 1]);
   ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Pres_500_spFreq_2'],'epsc');

fh=fh+1;
ind = [2:4:10];
for i=1:length(ind)
    figure(fh)
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
    
%     subplot(3,1,i)
    if i==length(ind)
        subplot(2,2,[3,4])
    else
        subplot(2,2,i)
    end
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');%,'Xlim',[0 0.6],'Ylim',[0 1]);
     ylim([0, 1]);   
     xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Pres_50_spFreq_10'],'epsc');


fh=fh+1;
ind = [4:4:12];
for i=1:length(ind)
    figure(fh)
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
%     subplot(3,1,i)
    if i==length(ind)
        subplot(2,2,[3,4])
    else
        subplot(2,2,i)
    end
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');%,'Xlim',[0 0.6],'Ylim',[0 1]);
    ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Pres_500_spFreq_10'],'epsc');

%% PROPORTION NO (FLIPPED PSYCHFIT beside original psychfit)
% fh=fh+1;
nBoots = 1;
fh=fh+1;
th_per_no = 1 - th_per;
ind = [1 3];
for i=1:length(ind)
    figure(fh)
    hold on
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
    misses = plotting(ind(i)).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
    
    subplot(2,2,2*i - 1)
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance', fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');
   ylim([0, 1]);
    xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
     
    subplot(2,2,2*i)
    psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
    ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc0_SpFreq_2_AlongsideFlip'],'epsc');

fh=fh+1;
ind = [2 4];
for i=1:length(ind)
    figure(fh)
    hold on
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
    misses = plotting(ind(i)).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
    
    subplot(2,2,2*i - 1)
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance', fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');
   ylim([0, 1]); 
   xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
     
    subplot(2,2,2*i)
    psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
   ylim([0, 1]);  
   xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc0_SpFreq_10_AlongsideFlip'],'epsc');

fh=fh+1;
ind = [5 7];
for i=1:length(ind)
    figure(fh)
    hold on
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
    misses = plotting(ind(i)).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
    
    subplot(2,2,2*i - 1)
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance', fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');
   ylim([0, 1]); 
   xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
     
    subplot(2,2,2*i)
    psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
   ylim([0, 1]);  
   xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc4_SpFreq_2_AlongsideFlip'],'epsc');

fh=fh+1;
ind = [6 8];
for i=1:length(ind)
    figure(fh)
    hold on
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
    misses = plotting(ind(i)).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
    
    subplot(2,2,2*i - 1)
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance', fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');
   ylim([0, 1]); 
   xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
     
    subplot(2,2,2*i)
    psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
   ylim([0, 1]);  
   xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc4_SpFreq_10_AlongsideFlip'],'epsc');

fh=fh+1;
ind = [9 11];
for i=1:length(ind)
    figure(fh)
    hold on
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
    misses = plotting(ind(i)).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
    
    subplot(2,2,2*i - 1)
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance', fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');
   ylim([0, 1]); 
   xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
     
    subplot(2,2,2*i)
    psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
   ylim([0, 1]);  
   xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc8_SpFreq_2_AlongsideFlip'],'epsc');

fh=fh+1;
ind = [10 12];
for i=1:length(ind)
    figure(fh)
    hold on
    lvl = plotting(ind(i)).contrast; % contrast
    hits = plotting(ind(i)).hits; % # correct
    misses = plotting(ind(i)).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(ind(i)).noOfStimDisp;
    fa_rate = plotting(ind(i)).falseAlarms(1)/plotting(ind(i)).lures(1);
    
    subplot(2,2,2*i - 1)
    psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance', fa_rate,'Lapses','Auto','Thresh',th_per,'DistType','Normal');
   ylim([0, 1]); 
   xlabel('Contrast')
     ylabel('Proportion Reported Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
     
    subplot(2,2,2*i)
    psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(ind(i)).ecc)...
        ,' ,Pres. - ',num2str(plotting(ind(i)).presTime),' ,Freq. - ',num2str(plotting(ind(i)).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
   ylim([0, 1]);  
   xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc8_SpFreq_10_AlongsideFlip'],'epsc');

%% PROPORTION NO (FLIPPED PSYCHFIT)
%draw eccentricity-wise, put the threshold value in the legend for now. Get
%the handle from the psyfit to adjust the color, legend, etc.
fh=fh+1;
nBoots = 1;

% fh=fh+1;
for i=1:4
    figure(fh)
    hold on
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
    misses = plotting(i).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    fa_rate = plotting(i).falseAlarms(1)/plotting(i).lures(1);
    subplot(2,2,i)
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    [thresh, par, threshB, parB]=psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
        ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
%     plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader
% %     plotting(i).thresh = round(mean(threshB),2);%Round because that's the resolution of the shader    
%     plotting(i).std = round(std(threshB),2);%Round because that's the resolution of the shader
%     ulvl = unique(lvl); % unique contrast levels
%     nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
% %     nlvl = plotting(i).noOfTrials;
%     ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
    ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc0_flip'],'epsc');

fh=fh+1;
for i=5:8
    figure(fh)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
    misses = plotting(i).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    fa_rate = plotting(i).falseAlarms(1)/plotting(i).lures(1);
    subplot(2,2,i-4)
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    [thresh, par, threshB, parB]=psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
        ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
    plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader
%     plotting(i).thresh = round(mean(threshB),2);%Round because that's the resolution of the shader    
    plotting(i).std = round(std(threshB),2);%Round because that's the resolution of the shader
    ulvl = unique(lvl); % unique contrast levels
    nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
%     nlvl = plotting(i).noOfTrials;
    ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
    ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc4_flip'],'epsc');

fh=fh+1;
for i=9:12
    figure(fh)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
    misses = plotting(i).misses;
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    fa_rate = plotting(i).falseAlarms(1)/plotting(i).lures(1);
    subplot(2,2,i-8)
%     [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
%         ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',1,'Extra');
    [thresh, par, threshB, parB]=psyfit(lvl,misses,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
        ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotHandle',fh,'Extra',...
       'Chance',fa_rate,'Lapses','Auto','Thresh',th_per_no,'Boots', nBoots, 'Flip','DistType','Normal');%'Xlim',[0 0.6],'Ylim',[0 1]);
    plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader
%     plotting(i).thresh = round(mean(threshB),2);%Round because that's the resolution of the shader    
    plotting(i).std = round(std(threshB),2);%Round because that's the resolution of the shader
    ulvl = unique(lvl); % unique contrast levels
    nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2); % number of trials at each contrast
%     nlvl = plotting(i).noOfTrials;
    ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./nlvl; % # hits at each contrast level
    ylim([0, 1]);
     xlabel('Contrast')
     ylabel('Proportion Reported NOT-Present')
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc8_flip'],'epsc');


%% D-prime

fh=fh+1;
for i=1:4
    figure(fh)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    subplot(2,2,i)
    plot(plotting(i).contrast, plotting(i).dprime,'LineWidth',2,'Color',[0 0.5 0.5])
    hold on
    %draws circles for displaying data points below
    for ii=1:length(plotting(i).contrast)
        plot( plotting(i).contrast(ii) , plotting(i).dprime(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7],'LineWidth',ceil(plotting(i).noOfStimDisp(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(plotting(i).contrast(ii) , plotting(i).dprime(ii), ...
            num2str(plotting(i).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
        hold on
    end
    %DISPLAY THRESHOLD AND TOTAL NUMBER OF TRIALS WITH STIMULUS DISPLAYED
    fontsize = 14;
    text( plotting(i).contrast(end) , 0,...
        strcat('threshold =',num2str(round(plotting(i).thresh*100)/100)),...
        'FontWeight','Bold','FontSize',fontsize,'Color',[.3 .3 .5]);
    text( plotting(i).contrast(end) , -0.5,...
        strcat('N = ',num2str(sum(plotting(i).noOfStimDisp))),...
        'FontSize',fontsize,'Color',[.3 .3 .5]);
    
    % Here we draw a line to mark the threshold
    [~,y2] = min(abs(plotting(i).contrast - plotting(i).thresh));
    %Vertical line
    line( [plotting(i).thresh plotting(i).thresh] , [min([plotting(i).dprime 0]) plotting(i).dprime(y2)], ...
        'color','r','LineStyle',':', 'LineWidth',2);
    %Horizontal line
    line( [0 plotting(i).thresh] , [plotting(i).dprime(y2) plotting(i).dprime(y2)],...
        'color','r','LineStyle',':', 'LineWidth',2);
    
     xlim([0, 0.6]);
     xlabel('Contrast')
     ylabel('D-Prime')
     title(['Ecc. - ',num2str(plotting(i).ecc),' ,Pres. - ',num2str(plotting(i).presTime),',Freq. - ',num2str(plotting(i).spFreq)])
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc0_Dprime'],'epsc');


fh = fh+1;
for i=5:8
    figure(fh)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    subplot(2,2,i-4)
    plot(plotting(i).contrast, plotting(i).dprime,'LineWidth',2,'Color',[0 0.5 0.5])
        hold on
    %draws circles for displaying data points below
    for ii=1:length(plotting(i).contrast)
        plot( plotting(i).contrast(ii) , plotting(i).dprime(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7],'LineWidth',ceil(plotting(i).noOfStimDisp(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(plotting(i).contrast(ii) , plotting(i).dprime(ii), ...
            num2str(plotting(i).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
        hold on
    end
    
    %DISPLAY THRESHOLD AND TOTAL NUMBER OF TRIALS WITH STIMULUS DISPLAYED
    fontsize = 14;
    text( plotting(i).contrast(end) , 0,...
        strcat('threshold =',num2str(round(plotting(i).thresh*100)/100)),...
        'FontWeight','Bold','FontSize',fontsize,'Color',[.3 .3 .5]);
    text( plotting(i).contrast(end) , -0.5,...
        strcat('N = ',num2str(sum(plotting(i).noOfStimDisp))),...
        'FontSize',fontsize,'Color',[.3 .3 .5]);
    
    % Here we draw a line to mark the threshold
    [~,y2] = min(abs(plotting(i).contrast - plotting(i).thresh));
    %Vertical line
    line( [plotting(i).thresh plotting(i).thresh] , [min([plotting(i).dprime 0]) plotting(i).dprime(y2)], ...
        'color','r','LineStyle',':', 'LineWidth',2);
    %Horizontal line
    line( [0 plotting(i).thresh] , [plotting(i).dprime(y2) plotting(i).dprime(y2)],...
        'color','r','LineStyle',':', 'LineWidth',2);
    
     xlim([0, 0.6]);
     xlabel('Contrast')
     ylabel('D-Prime')
     title(['Ecc. - ',num2str(plotting(i).ecc),' ,Pres. - ',num2str(plotting(i).presTime),',Freq. - ',num2str(plotting(i).spFreq)])
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc4_Dprime'],'epsc');

fh = fh+1;
for i=9:12
    figure(fh)
    lvl = plotting(i).contrast; % contrast
    hits = plotting(i).hits; % # correct
%     tr = plotting(i).noOfTrials;
    tr = plotting(i).noOfStimDisp;
    subplot(2,2,i-8)
    plot(plotting(i).contrast, plotting(i).dprime,'LineWidth',2,'Color',[0 0.5 0.5])
    hold on
    %draws circles for displaying data points below
    for ii=1:length(plotting(i).contrast)
        plot( plotting(i).contrast(ii) , plotting(i).dprime(ii),...
            'o','MarkerSize',20,'Color',...
            [.3 .8 .7],'LineWidth',ceil(plotting(i).noOfStimDisp(ii)/10));
        %SETS NUMBERS INSIDE DATA POINTS NK
        text(plotting(i).contrast(ii) , plotting(i).dprime(ii), ...
            num2str(plotting(i).noOfStimDisp(ii)),'Color',[1 .0 .0],...
            'HorizontalAlignment','Center');
        hold on
    end
    
   %DISPLAY THRESHOLD AND TOTAL NUMBER OF TRIALS WITH STIMULUS DISPLAYED
    fontsize = 14;
    text( plotting(i).contrast(end) , 0.5,...
        strcat('threshold =',num2str(round(plotting(i).thresh*100)/100)),...
        'FontWeight','Bold','FontSize',fontsize,'Color',[.3 .3 .5]);
    text( plotting(i).contrast(end) , 0,...
        strcat('N = ',num2str(sum(plotting(i).noOfStimDisp))),...
        'FontSize',fontsize,'Color',[.3 .3 .5]);
    
    % Here we draw a line to mark the threshold
    [~,y2] = min(abs(plotting(i).contrast - plotting(i).thresh));
    %Vertical line
    line( [plotting(i).thresh plotting(i).thresh] , [min([plotting(i).dprime 0]) plotting(i).dprime(y2)], ...
        'color','r','LineStyle',':', 'LineWidth',2);
    %Horizontal line
    line( [0 plotting(i).thresh] , [plotting(i).dprime(y2) plotting(i).dprime(y2)],...
        'color','r','LineStyle',':', 'LineWidth',2);
        
     xlim([0, 0.6]);
     xlabel('Contrast')
     ylabel('D-Prime')
     title(['Ecc. - ',num2str(plotting(i).ecc),' ,Pres. - ',num2str(plotting(i).presTime),',Freq. - ',num2str(plotting(i).spFreq)])
     set(gca,'box','off','tickdir','out','FontSize',20)
end
% set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
saveas(figure(fh),[sub_fig_path '\Ecc8_Dprime'],'epsc');

% %% Sensitivity/Threshold summary plots
% %Pres time vs sensitivity (2 cpd)
% ecc_levels = [0 4 8];
% fh=fh+1;
% figure(fh)
% alpha = [];
% % for c=1:length(freq_levels)
% %     for p = 1:length(pres_levels)
%         for i=1:length(ecc_levels)
%             In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
%             In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
%             t1 = plotting(In1).thresh;
%             t2 = plotting(In2).thresh;
%             std1 = plotting(In1).std;
%             std2 = plotting(In2).std;
%             stddev = [1/std1 1/std2];
%             if t1>0.6
%                 t1 = 0.6;
%             elseif t1<0
%                 t1 = 0.01;
%             end
%             if t2>0.6
%                 t2 = 0.6;
%             elseif t2<0
%                 t2 = 0.01;
%             end
%             sens = [1/t1 1/t2];
%             alpha = [alpha (sens(2) - sens(1))/(450)];
%             formatSpec = "(sens/ms)";
% %             if sens>0.5
% %                 sens = 0.5;
% %             elseif sens<0
% %                 sens = 0;
% %             end
%             if ecc_levels(i)==0
%                 plot([50 500],sens,'-s','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
% %                 errorbar([50 500],sens,stddev,'-s','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
%                 txt = ['\angle' alpha(i) formatSpec];
%                 %                 txt = ['\angle',alpha(i)];
%                 text(50,sens(1),txt,'FontSize',20)
%                 hold on
%             elseif ecc_levels(i)==4
%                 plot([50 500],sens,'-v','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
% %                 errorbar([50 500],sens,stddev,'-v','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
%                 txt = ['\angle' alpha(i) formatSpec];
%                 %                 txt = ['\angle',alpha(i)];
%                 text(50,sens(1),txt,'FontSize',20)
%                 hold on
%             elseif ecc_levels(i)==8
%                 plot([50 500],sens,'-o','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
% %                 errorbar([50 500],sens,stddev,'-o','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
%                 txt = ['\angle' alpha(i) formatSpec];
%                 %                 txt = ['\angle',alpha(i)];
%                 text(50,sens(1),txt,'FontSize',20)
%                 hold on
% %             elseif ecc_levels(i)==0 && freq_levels(c)==10
% %                 plot([50 500],th,'--s','linewidth', 2,'Color','k','MarkerEdgeColor','k',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==4 && freq_levels(c)==10
% %                 plot([50 500],th,'--v','linewidth', 2,'Color','r','MarkerEdgeColor','r',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==8 && freq_levels(c)==10
% %                 plot([50 500],th,'--o','linewidth', 2,'Color','b','MarkerEdgeColor','b',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
%             end
%         end
% %     end
%     
% % end
% hold off
% set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% legend
% xticks([50 500])
% xlim([-50 600])
% xlabel('Presentation time')
% ylabel('Sensitivity')
% title('Low SpFreq')
% saveas(figure(fh),[sub_fig_path '\Sens_Pres_2'],'epsc');
% 
% %Pres time vs sensitivity (10 cpd)
% 
% alpha = [];
% fh=fh+1;
% figure(fh)
% % y_tick=[];
% % for c=1:length(freq_levels)
% %     for p = 1:length(pres_levels)
%         for i=1:length(ecc_levels)
%             In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
%             In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
%             t1 = plotting(In1).thresh;
%             t2 = plotting(In2).thresh;
%             std1 = plotting(In1).std;
%             std2 = plotting(In2).std;
%             stddev = [1/std1 1/std2];
%             if t1>0.6
%                 t1 = 0.6;
%             elseif t1<0
%                 t1 = 0.01;
%             end
%             if t2>0.6
%                 t2 = 0.6;
%             elseif t2<0
%                 t2 = 0.01;
%             end
%             sens = [1/t1 1/t2];
% %             alpha = [alpha atand((sens(2) - sens(1))/(450))];            
%             alpha = [alpha (sens(2) - sens(1))/(450)];
%             formatSpec = "(sens/ms)";
% %             if ecc_levels(i)==0 && freq_levels(c)==2
% %                 plot([50 500],th,'-s','linewidth', 2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==4 && freq_levels(c)==2
% %                 plot([50 500],th,'-v','linewidth', 2,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==8 && freq_levels(c)==2
% %                 plot([50 500],th,'-o','linewidth', 2,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
%             if ecc_levels(i)==0
%                 plot([50 500],sens,'-s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
% %                 errorbar([50 500],sens,stddev,'-s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 %                 txt = sprintf(formatSpec,alpha(i));
%                 txt = ['\angle' alpha(i) formatSpec];
% %                 txt = ['\angle',alpha(i)];
%                 text(50,sens(1),txt,'FontSize',20)
%                 hold on
%             elseif ecc_levels(i)==4
% %                 errorbar([50 500],sens,stddev,'-v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 plot([50 500],sens,'-v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 txt = ['\angle' alpha(i) formatSpec];
%                 text(50,sens(1),txt,'FontSize',20)
%                 hold on
%             elseif ecc_levels(i)==8
% %                 errorbar([50 500],sens,stddev,'-o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 plot([50 500],sens,'-o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 txt = ['\angle' alpha(i) formatSpec];
%                 text(50,sens(1),txt,'FontSize',20)
%                 hold on
%             end
% %             y_tick = [y_tick th];
%         end
% %     end
% %     
% % end
% hold off
% set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% legend
% xticks([50 500])
% xlim([-50 600])
% % ylim([0 0.5])
% % yticks(sort(y_tick))
% xlabel('Presentation time')
% ylabel('Sensitivity')
% title('High SpFreq')
% saveas(figure(fh),[sub_fig_path '\Sens_Pres_10'],'epsc');
% %%
% 
% %Pres time vs threshold (2 cpd)
% fh=fh+1;
% figure(fh)
% % for c=1:length(freq_levels)
% %     for p = 1:length(pres_levels)
%         for i=1:length(ecc_levels)
%             In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
%             In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
%             t1 = plotting(In1).thresh;
%             t2 = plotting(In2).thresh;
%             std1 = plotting(In1).std;
%             std2 = plotting(In2).std;
%             stddev = [std1 std2];
%             if t1>0.6
%                 t1 = 0.6;
%             elseif t1<0
%                 t1 = 0.01;
%             end
%             if t2>0.6
%                 t2 = 0.6;
%             elseif t2<0
%                 t2 = 0.01;
%             end
%             thresh = [t1 t2];
%             if ecc_levels(i)==0
%                 errorbar([50 500],thresh,stddev,'--s','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
%                 hold on
%             elseif ecc_levels(i)==4
%                 errorbar([50 500],thresh,stddev,'--v','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
%                 hold on
%             elseif ecc_levels(i)==8
%                 errorbar([50 500],thresh,stddev,'--o','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
%                 hold on
% %             elseif ecc_levels(i)==0 && freq_levels(c)==10
% %                 plot([50 500],th,'--s','linewidth', 2,'Color','k','MarkerEdgeColor','k',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==4 && freq_levels(c)==10
% %                 plot([50 500],th,'--v','linewidth', 2,'Color','r','MarkerEdgeColor','r',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==8 && freq_levels(c)==10
% %                 plot([50 500],th,'--o','linewidth', 2,'Color','b','MarkerEdgeColor','b',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
%             end
%         end
% %     end
%     
% % end
% hold off
% set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% legend
% xticks([50 500])
% xlim([-50 600])
% ylim([0 0.6])
% xlabel('Presentation time')
% ylabel('Threshold')
% title('Low SpFreq')
% saveas(figure(fh),[sub_fig_path '\Thresh_Pres_2'],'epsc');
% 
% %Pres time vs threshold (10 cpd)
% fh=fh+1;
% figure(fh)
% % y_tick=[];
% % for c=1:length(freq_levels)
% %     for p = 1:length(pres_levels)
%         for i=1:length(ecc_levels)
%             In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
%             In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
%             t1 = plotting(In1).thresh;
%             t2 = plotting(In2).thresh;
%             std1 = plotting(In1).std;
%             std2 = plotting(In2).std;
%             stddev = [std1 std2];
%             if t1>0.6
%                 t1 = 0.6;
%             elseif t1<0
%                 t1 = 0.01;
%             end
%             if t2>0.6
%                 t2 = 0.6;
%             elseif t2<0
%                 t2 = 0.01;
%             end
%             thresh = [t1 t2];
%             
% %             if ecc_levels(i)==0 && freq_levels(c)==2
% %                 plot([50 500],th,'-s','linewidth', 2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==4 && freq_levels(c)==2
% %                 plot([50 500],th,'-v','linewidth', 2,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
% %             elseif ecc_levels(i)==8 && freq_levels(c)==2
% %                 plot([50 500],th,'-o','linewidth', 2,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
% %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
% %                 hold on
%             if ecc_levels(i)==0
%                 errorbar([50 500],thresh,stddev,'--s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 hold on
%             elseif ecc_levels(i)==4
%                 errorbar([50 500],thresh,stddev,'--v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 hold on
%             elseif ecc_levels(i)==8
%                 errorbar([50 500],thresh,stddev,'--o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
%                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
%                 hold on
%             end
% %             y_tick = [y_tick th];
%         end
% %     end
% %     
% % end
% hold off
% set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% legend
% xticks([50 500])
% xlim([-50 600])
% ylim([0 0.6])
% % yticks(sort(y_tick))
% xlabel('Presentation time')
% ylabel('Threshold')
% title('High SpFreq')
% saveas(figure(fh),[sub_fig_path '\Thresh_Pres_10'],'epsc');
% %% Radial Colorbar plot of thresholds
% 
% imsize_x = 1080;
% imsize_y = 1920;
% pixelAngle = 1.06; % arcmin per pixel, was 0.8 initially
% width = 1;
% fh=fh+1;
% fig_count =fh;
% for p = 1:length(pres_levels) %different presentation times
%     for c=1:length(freq_levels)
%         grat = zeros(imsize_x,imsize_y);        
%         for ii = 1:imsize_x
%             for jj = 1:imsize_y
%                 r = sqrt ((ii - imsize_x/2)^2 + (jj - imsize_y/2)^2)* pixelAngle / 60;
%                 if ((r >= ecc_levels(1) - width/2) && (r <= ecc_levels(1) + width/2))
%                     In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(1) & [plotting.spFreq]==freq_levels(c)); %indices
%                     scale = plotting(In1).thresh;
%                 elseif ((r >= ecc_levels(2) - width/2) && (r <= ecc_levels(2) + width/2))
%                     In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(2) & [plotting.spFreq]==freq_levels(c)); %indices
%                     scale = plotting(In1).thresh;
%                 elseif ((r >= ecc_levels(3) - width/2) && (r <= ecc_levels(3) + width/2))
%                     In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(3) & [plotting.spFreq]==freq_levels(c)); %indices
%                     scale = plotting(In1).thresh;
%                 else
%                     scale = 0;
%                 end
%                 grat(ii, jj) = scale;
%             end
%         end
%         figure(fig_count)
%         imagesc(grat);
%         set (gca, 'CLim',[0 0.5])
%         colorbar
%         title(['Pres. - ',num2str(pres_levels(p)),' ,Freq. - ',num2str(freq_levels(c))])
%         saveas(figure(fig_count),[sub_fig_path '\Color_' ['Pres_',num2str(pres_levels(p)),'_Freq_',num2str(freq_levels(c))]],'epsc');
%         fig_count =fig_count+1;
%     end
% end
% 
% 
