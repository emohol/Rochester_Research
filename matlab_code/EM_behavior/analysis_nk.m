%EM_analysis: Amplitude, angle, duration distributions

path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
% subject = 'A013';%'Nikunj'
% subject = 'Nikunj';
subject = 'A092';
% subject = 'A036';
pathtodata = fullfile(path,subject);

bin = 0; %set to 1 if binning
read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

%%
fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures';
sub_fig_path = fullfile(fig_path,subject,'EM');
if (~isfolder(sub_fig_path))
    mkdir(sub_fig_path)
end

[valid,counter] = countingTrialsNK(ppt);

trials = struct;
MS = struct;
S = struct;
P_S = struct; %Post saccade saccades
count_MS = 1;
count_S = 1;
count_PS = 1;

for i=1:length(ppt)
    if valid.validTrials(i)
%     if valid.drift(i)
        timeOn = round(ppt{i}.saccOn);%Exposure period
        stimOff = round(ppt{i}.stimOff);
        [Answer, Intersected] = isIntersectedIn(timeOn,stimOff-timeOn,ppt{i}.microsaccades);
        if Answer
            MS(count_MS).amp =  ppt{i}.microsaccades.amplitude(Intersected);
            MS(count_MS).angle =  ppt{i}.microsaccades.angle(Intersected);
            count_MS = count_MS+1;
        end
        
        timeOn = round(ppt{i}.saccOn);%Saccade period
        timeOff = round(ppt{i}.saccOff);
        [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,ppt{i}.saccades);
        if Answer
            S(count_S).amp =  ppt{i}.saccades.amplitude(Intersected);
            S(count_S).angle =  ppt{i}.saccades.angle(Intersected);
            S(count_S).dur =  timeOff-timeOn;
            S(count_S).onSpeed = ppt{i}.velocity(timeOn);
%             S(count_S).vel = ppt{i}.velocity(timeOn:timeOff);            
            S(count_S).dur_3 = sum(ppt{i}.velocity(timeOn:timeOff) < 180);
            S(count_S).landDist = sqrt(power((ppt{i}.x.position(timeOff)+ ppt{i}.xoffset * ppt{i}.pixelAngle),2)...
                + power((ppt{i}.y.position(timeOff)+ ppt{i}.yoffset * ppt{i}.pixelAngle),2));
            count_S = count_S+1;
        end
        
        [~,emat] = (min(abs(ppt{i}.drifts.start - timeOff)));
        if ~isempty(emat)
            EMAT_driftOn = ppt{i}.drifts.start(emat);
            EMAT_driftOff = EMAT_driftOn + ppt{i}.drifts.duration(emat);
        end
        timeOn = round(ppt{i}.saccOff);%Post-Saccade period
        timeOff = round(ppt{i}.stimOff);
        [Answer, Intersected] = isIntersectedIn(EMAT_driftOn,stimOff-EMAT_driftOn,ppt{i}.saccades);
%         [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,ppt{i}.saccades);
        if Answer
            P_S(count_PS).amp =  ppt{i}.saccades.amplitude(Intersected);
            P_S(count_PS).angle =  ppt{i}.saccades.angle(Intersected);
            count_PS = count_PS+1;
        end
%          for( k = 1 : size(ppt{i}.microsaccades.start,2) )
%                 if(ppt{i}.microsaccades.start(k) >= ppt{i}.saccOff &&...
%                        ppt{i}.microsaccades.start(k) < ppt{i}.stimOff)
%                     MS(count).amp =  ppt{i}.microsaccades.amplitude(k);
%                     MS(count).angle =  ppt{i}.microsaccades.angle(k);
%                     count = count+1;
%                 end
%          end
%         [minValue,closestIndex]=min(abs(ppt{i}.saccOff - ppt{i}.microsaccades.start));
%         if minValue~=0
%             MS(count).amp =  ppt{i}.microsaccades.amplitude(closestIndex);
%             MS(count).angle =  ppt{i}.microsaccades.angle(closestIndex);
%             count = count+1;
%         end
    end
end
allMSAmps = cat(2, MS.amp);
ntotal = length(allMSAmps);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(0, 60, nbins);
% c = linspace(min(allMSAmps), max(allMSAmps), nbins);

m1 = nanmean(allMSAmps(:));
s1 = nanstd(allMSAmps(:));% / sqrt(sum(~isnan(allMSAmps(:))));
n1 = hist(allMSAmps(:), c);
n1 = n1 / sum(n1);

% out.(sprintf('amp_data')) = {allMSAmps(:)};
% out.(sprintf('amp_hist')) = {n1};
% out.(sprintf('amp_hist_bins')) = {c};
% out.(sprintf('amp_MS')) = {[m1, s1]};
fh = 1;
figure(fh)
hist(allMSAmps(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('amplitude (arcmin)');
ylabel('Micro-sacc prob');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(m1,max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\MS-Amp','epsc');
saveas(figure(fh),[sub_fig_path '\MS-Amp'],'epsc');

fh=fh+1;
allMSAng = cat(2, MS.angle);
% ntotal = length(allMSAng);
% nbins = max(10, min(60, round(ntotal/3)));
c = linspace(0, 2*pi, 30);
n1 = hist(allMSAng(:), c);
n1 = n1 / sum(n1);%Convert to prob
figure(fh)
h=polar([c, c(1)], [n1, n1(1)]); % c(1) and n1(1) to complete the circle
set(h, 'linewidth', 2);
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Micro-sacc direction prob');
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\MS-Ang','epsc');
saveas(figure(fh),[sub_fig_path '\MS-Ang'],'epsc');

fh=fh+1;
allSacAmps = cat(2, S.amp);
ntotal = length(allSacAmps);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacAmps), max(allSacAmps), nbins);
m1 = nanmean(allSacAmps(:));
s1 = nanstd(allSacAmps(:));% / sqrt(sum(~isnan(allSacAmps(:))));
n1 = hist(allSacAmps(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacAmps(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('amplitude (arcmin)');
ylabel('Sacc prob');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(m1,max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-Amp','epsc');
saveas(figure(fh),[sub_fig_path '\Sac-Amp'],'epsc');

fh=fh+1;
allSacDurs = cat(2, S.dur);
ntotal = length(allSacDurs);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacDurs), max(allSacDurs), nbins);
m1 = nanmean(allSacDurs(:));
s1 = nanstd(allSacDurs(:));% / sqrt(sum(~isnan(allSacDurs(:))));
n1 = hist(allSacDurs(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacDurs(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Duration (ms)');
ylabel('Prob');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(m1,max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-Dur','epsc');
saveas(figure(fh),[sub_fig_path '\Sac-Dur'],'epsc');

fh=fh+1;
allSacOnVel = cat(2, S.onSpeed);
allSacOnVel = allSacOnVel./60;
ntotal = length(allSacOnVel);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacOnVel), max(allSacOnVel), nbins);
m1 = nanmean(allSacOnVel(:));
s1 = nanstd(allSacOnVel(:));% / sqrt(sum(~isnan(allSacOnVel(:))));
n1 = hist(allSacOnVel(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacOnVel(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Onset Speed (deg/s)');
ylabel('Prob');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-OnVel','epsc');
saveas(figure(fh),[sub_fig_path '\Sac-OnVel'],'epsc');

fh=fh+1;
allSac3dur = cat(2, S.dur_3);
ntotal = length(allSac3dur);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSac3dur), max(allSac3dur), nbins);
m1 = nanmean(allSac3dur(:));
s1 = nanstd(allSac3dur(:)) ;%/ sqrt(sum(~isnan(allSac3dur(:))));
n1 = hist(allSac3dur(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSac3dur(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Duration with speed < 180 arcmin/s (ms)');
ylabel('Prob');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-Dur-3','epsc');
saveas(figure(fh),[sub_fig_path '\Sac-Dur-3'],'epsc');
% dur_3 = sum(ppt{end}.velocity(timeOn:timeOff) < 180);

%Currently saccade landing distance is analysed by condition in
%'analysis_by_condition.m'
% fh=fh+1;
% allSacLandDis = cat(2, S.landDist);
% ntotal = length(allSacLandDis);
% nbins = max(10, min(60, round(ntotal/3)));
% c = linspace(min(allSacLandDis), max(allSacLandDis), nbins);
% m1 = nanmean(allSacLandDis(:));
% s1 = nanstd(allSacLandDis(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
% n1 = hist(allSacLandDis(:), c);
% n1 = n1 / sum(n1);
% figure(fh)
% hist(allSacLandDis(:), c);
% % plot(c, n1,'linewidth', 2);
% vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
% set(gca,'box','off','tickdir','out','FontSize',20)
% % Enlarge figure to full screen.
% % set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% xlabel('Landing distance (arcmin)');
% ylabel('Prob');
% formatSpec = "Mean - %1.1f, Std - %1.1f";
% txt = sprintf(formatSpec,m1,s1);
% text(double(m1),max(n1),txt,'FontSize',20);
% % saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-Land','epsc');
% saveas(figure(fh),[sub_fig_path '\Sac-Land'],'epsc');

%Post saccade saccades
allPSAmps = cat(2, P_S.amp);
ntotal = length(allPSAmps);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(0, 60, nbins);
m1 = nanmean(allPSAmps(:));
s1 = nanstd(allPSAmps(:));% / sqrt(sum(~isnan(allMSAmps(:))));
n1 = hist(allPSAmps(:), c);
n1 = n1 / sum(n1);
fh=fh+1;
figure(fh)
% hist(allPSAmps(:), c);
hist(allPSAmps(:));
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('amplitude (arcmin)');
ylabel('Post-sacc saccades prob');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(m1,max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\PS-Amp','epsc');
saveas(figure(fh),[sub_fig_path '\PS-Amp'],'epsc');

fh=fh+1;
allPSAng = cat(2, P_S.angle);
% ntotal = length(allMSAng);
% nbins = max(10, min(60, round(ntotal/3)));
c = linspace(0, 2*pi, 30);
n1 = hist(allPSAng(:), c);
n1 = n1 / sum(n1);%Convert to prob
figure(fh)
h=polar([c, c(1)], [n1, n1(1)]); % c(1) and n1(1) to complete the circle
set(h, 'linewidth', 2);
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Post-sacc saccades direction prob');
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\PS-Ang','epsc');
saveas(figure(fh),[sub_fig_path '\PS-Ang'],'epsc');
%% Drift


      

