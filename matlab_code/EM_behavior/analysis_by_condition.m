%EM_analysis: Landing distance by eccentricity
path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
% subject = 'A013';
% subject = 'Nikunj';
% subject = 'A092';
subject = 'A036';
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



% [valid,counter] = countingTrialsNK(ppt);
[valid,counter] = countingTrialsNK_noLand(ppt);

trials = struct;
% MS = struct;
% S = struct;
% P_S = struct; %Post saccade saccades
% count_MS = 1;
% count_S = 1;
% count_PS = 1;

E0= struct;
E4 = struct;
E8 = struct;
count_E0 = 1;
count_E4 = 1;
count_E8 = 1;
for i=1:length(ppt)
%     if valid.validTrials(i)
    if valid.drift(i)
        timeOn = round(ppt{i}.saccOn);%Saccade period
        timeOff = round(ppt{i}.saccOff);
        [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,ppt{i}.saccades);
        if Answer
            if (ppt{i}.eccentricity==0)
                E0(count_E0).x = ppt{i}.x.position(timeOff)+ ppt{i}.xoffset * ppt{i}.pixelAngle;
                E0(count_E0).y = ppt{i}.y.position(timeOff)+ ppt{i}.yoffset * ppt{i}.pixelAngle;
                E0(count_E0).landDist = sqrt(power((ppt{i}.x.position(timeOff)+ ppt{i}.xoffset * ppt{i}.pixelAngle),2)...
                    + power((ppt{i}.y.position(timeOff)+ ppt{i}.yoffset * ppt{i}.pixelAngle),2));
                E0(count_E0).angle = atan2(E0(count_E0).y,E0(count_E0).x);
                if E0(count_E0).angle<0
                    E0(count_E0).angle = E0(count_E0).angle + 360;
                end
                count_E0 = count_E0+1;
                
            elseif (ppt{i}.eccentricity==4)
                E4(count_E4).x = ppt{i}.x.position(timeOff)+ ppt{i}.xoffset * ppt{i}.pixelAngle;
                E4(count_E4).y = ppt{i}.y.position(timeOff)+ ppt{i}.yoffset * ppt{i}.pixelAngle;
                E4(count_E4).landDist = sqrt(power((ppt{i}.x.position(timeOff)+ ppt{i}.xoffset * ppt{i}.pixelAngle),2)...
                    + power((ppt{i}.y.position(timeOff)+ ppt{i}.yoffset * ppt{i}.pixelAngle),2));
                E4(count_E4).angle = atan2(E4(count_E4).y,E4(count_E4).x);
                if E4(count_E4).angle<0
                    E4(count_E4).angle = E4(count_E4).angle + 360;
                end
                count_E4 = count_E4+1;
                
            elseif (ppt{i}.eccentricity==8)
                E8(count_E8).x = ppt{i}.x.position(timeOff)+ ppt{i}.xoffset * ppt{i}.pixelAngle;
                E8(count_E8).y = ppt{i}.y.position(timeOff)+ ppt{i}.yoffset * ppt{i}.pixelAngle;
                E8(count_E8).landDist = sqrt(power((ppt{i}.x.position(timeOff)+ ppt{i}.xoffset * ppt{i}.pixelAngle),2)...
                    + power((ppt{i}.y.position(timeOff)+ ppt{i}.yoffset * ppt{i}.pixelAngle),2));
                E8(count_E8).angle = atan2(E8(count_E8).y,E8(count_E8).x);
                if E8(count_E8).angle<0
                    E8(count_E8).angle = E8(count_E8).angle + 360;
                end
                count_E8 = count_E8+1;
            end
        end
    end
end
%%
fh=1;
allSacLandDis = cat(2, E0.landDist);
ntotal = length(allSacLandDis);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis), max(allSacLandDis), nbins);
m1 = nanmean(allSacLandDis(:));
s1 = nanstd(allSacLandDis(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance (arcmin)');
ylabel('Prob');
title('Ecc 0');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-Land','epsc');
saveas(figure(fh),[sub_fig_path '\Ecc_0_Land'],'epsc');
fh=fh+1;

allSacLandDis_x = cat(2, E0.x);
ntotal = length(allSacLandDis_x);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis_x), max(allSacLandDis_x), nbins);
m1 = nanmean(allSacLandDis_x(:));
s1 = nanstd(allSacLandDis_x(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis_x(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis_x(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance X(arcmin)');
ylabel('Prob');
title('Ecc 0');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
saveas(figure(fh),[sub_fig_path '\Ecc_0_Land_x'],'epsc');
fh=fh+1;


allSacLandDis_y = cat(2, E0.y);
ntotal = length(allSacLandDis_y);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis_y), max(allSacLandDis_y), nbins);
m1 = nanmean(allSacLandDis_y(:));
s1 = nanstd(allSacLandDis_y(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis_y(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis_y(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance Y(arcmin)');
ylabel('Prob');
title('Ecc 0');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
saveas(figure(fh),[sub_fig_path '\Ecc_0_Land_y'],'epsc');
fh=fh+1;

allMSAng = cat(2, E0.angle);
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
xlabel('Direction prob');
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\MS-Ang','epsc');
saveas(figure(fh),[sub_fig_path '\Ecc_0_Land_Polar'],'epsc');
fh=fh+1;

allSacLandDis = [];
allSacLandDis = cat(2, E4.landDist);
ntotal = length(allSacLandDis);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis), max(allSacLandDis), nbins);
m1 = nanmean(allSacLandDis(:));
s1 = nanstd(allSacLandDis(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance (arcmin)');
ylabel('Prob');
title('Ecc 4');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-Land','epsc');
saveas(figure(fh),[sub_fig_path '\Ecc_4_Land'],'epsc');
fh=fh+1;

allSacLandDis_x = [];
allSacLandDis_x = cat(2, E4.x);
ntotal = length(allSacLandDis_x);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis_x), max(allSacLandDis_x), nbins);
m1 = nanmean(allSacLandDis_x(:));
s1 = nanstd(allSacLandDis_x(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis_x(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis_x(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance X(arcmin)');
ylabel('Prob');
title('Ecc 4');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
saveas(figure(fh),[sub_fig_path '\Ecc_4_Land_x'],'epsc');
fh=fh+1;

allSacLandDis_y = [];
allSacLandDis_y = cat(2, E4.y);
ntotal = length(allSacLandDis_y);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis_y), max(allSacLandDis_y), nbins);
m1 = nanmean(allSacLandDis_y(:));
s1 = nanstd(allSacLandDis_y(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis_y(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis_y(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance Y(arcmin)');
ylabel('Prob');
title('Ecc 4');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
saveas(figure(fh),[sub_fig_path '\Ecc_4_Land_y'],'epsc');
fh=fh+1;

allMSAng = cat(2, E4.angle);
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
xlabel('Direction prob');
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\MS-Ang','epsc');
saveas(figure(fh),[sub_fig_path '\Ecc_4_Land_Polar'],'epsc');
fh=fh+1;

allSacLandDis = [];
allSacLandDis = cat(2, E8.landDist);
ntotal = length(allSacLandDis);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis), max(allSacLandDis), nbins);
m1 = nanmean(allSacLandDis(:));
s1 = nanstd(allSacLandDis(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance (arcmin)');
ylabel('Prob');
title('Ecc 8');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\Sac-Land','epsc');
saveas(figure(fh),[sub_fig_path '\Ecc_8_Land'],'epsc');
fh=fh+1;


allSacLandDis_x = cat(2, E8.x);
ntotal = length(allSacLandDis_x);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis_x), max(allSacLandDis_x), nbins);
m1 = nanmean(allSacLandDis_x(:));
s1 = nanstd(allSacLandDis_x(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis_x(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis_x(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance X(arcmin)');
ylabel('Prob');
title('Ecc 8');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
saveas(figure(fh),[sub_fig_path '\Ecc_8_Land_x'],'epsc');
fh=fh+1;


allSacLandDis_y = cat(2, E8.y);
ntotal = length(allSacLandDis_y);
nbins = max(10, min(60, round(ntotal/3)));
c = linspace(min(allSacLandDis_y), max(allSacLandDis_y), nbins);
m1 = nanmean(allSacLandDis_y(:));
s1 = nanstd(allSacLandDis_y(:));% / sqrt(sum(~isnan(allSacLandDis(:))));
n1 = hist(allSacLandDis_y(:), c);
n1 = n1 / sum(n1);
figure(fh)
hist(allSacLandDis_y(:), c);
% plot(c, n1,'linewidth', 2);
vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
set(gca,'box','off','tickdir','out','FontSize',20)
% Enlarge figure to full screen.
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('Landing distance Y(arcmin)');
ylabel('Prob');
title('Ecc 8');
formatSpec = "Mean - %1.1f, Std - %1.1f";
txt = sprintf(formatSpec,m1,s1);
text(double(m1),max(n1),txt,'FontSize',20);
saveas(figure(fh),[sub_fig_path '\Ecc_8_Land_y'],'epsc');
fh=fh+1;

allMSAng = cat(2, E8.angle);
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
xlabel('Direction prob');
% saveas(figure(fh),'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\EM\MS-Ang','epsc');
saveas(figure(fh),[sub_fig_path '\Ecc_8_Land_Polar'],'epsc');
fh=fh+1;