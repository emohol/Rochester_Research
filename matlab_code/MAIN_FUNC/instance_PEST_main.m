function [figure(1)] = instance_PEST_main(ppt)
trials = struct;
for i=1:noFiles
%     trials.contrast(i) = ppt{i,1}.contrast;
%     trials.resp(i) = ppt{i,1}.resp;
%     trials.present(i) = ppt{i,1}.present;
    trials(i).contrast = ppt{i,1}.contrast;
    trials(i).resp = ppt{i,1}.resp;
    trials(i).present = ppt{i,1}.present;
    trials(i).respTime = ppt{i,1}.responseTime;
    trials(i).ecc = ppt{i,1}.eccentricity;
    trials(i).presTime = ppt{i,1}.presTime;
    trials(i).spFreq = ppt{i,1}.spatialFreq;
    
    trials(i).pixelAngle = ppt{i,1}.pixelAngle;
    trials(i).xoffset = ppt{i,1}.xoffset;
    trials(i).yoffset = ppt{i,1}.yoffset;
    trials(i).cuePrf = ppt{i,1}.cuePrf;
    trials(i).fixOn = ppt{i,1}.fixOn;
    trials(i).presTime = ppt{i,1}.presTime;
    
    trials(i).cueCtr = ppt{i,1}.cueCtr;
    trials(i).saccOn = ppt{i,1}.saccOn;
    trials(i).saccOff = ppt{i,1}.saccOff;
    trials(i).flashOn = ppt{i,1}.flashOn;
    trials(i).rampOff = ppt{i,1}.rampOff;
    trials(i).stimOff = ppt{i,1}.stimOff;
    trials(i).quit = ppt{i,1}.quit; 
    
    %Correct for erroenous PEST levels
    if trials(i).contrast == 0 && trials(i).present == 1
        trials(i).present = 0;
    end
    
    if trials(i).contrast > 0.6
        trials(i).contrast = 0.6;
    end
end

%round off to 'bin' contrast values
% for i=1:noFiles
%     trials(i).contrast = round([trials(i).contrast],1);
% end

% boxplot(x,{H,M})
% grpstats(x,{H,M})
% [u g]=grpstats(x,{H,M},{'mean','gname'})
% [u g n]=grpstats(x,{H,M},{'mean','gname','numel'})

ecc_levels = unique([trials.ecc]);
% con_lev= unique(contrast(In2));
contrast = [trials.contrast];
pres_Time = [trials.presTime];
sp_Freq = [trials.spFreq];
% x=sort([trials.contrast]);
T = table;
count = 1;
for i=1:length(ecc_levels)
    plotting = struct([]);
% for i=length(ecc_levels)
% for i=1:1:2
% for i=2:1:3  
    In1=([trials.ecc]==ecc_levels(i)); %indices
%     con_lev= unique(contrast(In1));    
%     In2=([trials.presTime]==pres_levels(i)); %indices
    pres_levels= unique(pres_Time(In1));
    [~,u] = size(pres_levels);
    for p = 1:u %different presentation times
        plotting(p).presTime = pres_levels(p);
        In2=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)); %indices
        spFreq_lev = unique(sp_Freq(In2));
%         con_lev= unique(contrast(In2));
        plotting(p).ecc = ecc_levels(i);
%         plotting(p) = struct('ecc',ecc_levels(i));
%         plotting(p).contrast = con_lev;
        
        [~,t] = size(spFreq_lev);
        for j = 1:t
            plotting(p).spFreq = spFreq_lev(j);
%             In3=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1); %indices
%             con_lev= contrast(In3);
%             x = find(In3);
            
            In3_resp1=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1 &[trials.resp]==1); %indices
            con_lev_resp1 = contrast(In3_resp1);
            x_resp1 = find(In3_resp1); %Convert Logical into Sequential Index Vector 
            
%             In3_resp0=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1 &[trials.resp]==0); %indices
            In3_resp0=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1); %indices
            con_lev_resp0 = contrast(In3_resp0);
            x_resp0 = find(In3_resp0); %Convert Logical into Sequential Index Vector 
            last_con=trials(x_resp0(end)).contrast;
            figure(1);
%             for sc=1:length(x_resp1)
%                 if trials(x_resp1(sc)).resp == 1%check if response is 1                    
%                     scatter(x_resp1(sc),con_lev_resp1(sc),'d','filled')
% %                     hold on
%                 else %response is 2
%                     scatter(x_resp1(sc),con_lev_resp1(sc),'d')
% %                     hold on
%                 end
%                 hold on
%             end
            if (spFreq_lev(j) == 2 && ecc_levels(i)== 0)
                han(count) = plot(x_resp0,con_lev_resp0,':d','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1))...
                    ,',Amp. - ',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (spFreq_lev(j) == 2 && ecc_levels(i)== 4)
                han(count) = plot(x_resp0,con_lev_resp0,'--d','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (spFreq_lev(j) == 2 && ecc_levels(i)== 8)
                han(count) = plot(x_resp0,con_lev_resp0,'-d','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (spFreq_lev(j) == 10 && ecc_levels(i)== 0)
                han(count) = plot(x_resp0,con_lev_resp0,':s','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))...
                   ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1))...
                   ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'s','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color) 
            elseif (spFreq_lev(j) == 10 && ecc_levels(i)== 4)
                han(count) = plot(x_resp0,con_lev_resp0,'--s','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'s','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (spFreq_lev(j) == 10 && ecc_levels(i)== 8)
                han(count) = plot(x_resp0,con_lev_resp0,'-s','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==spFreq_lev(j) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'s','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            end
            count = count+1;
%             plot(x_resp1,con_lev_resp1,'d-','MarkerFaceColor','black','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))])
%             hold on
%             plot(x_resp0,con_lev_resp0,'d-','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(spFreq_lev(j))])
        end        
    end    
%     plotting = struct([]);
end

In4 = ([trials.present]==0);
x = find(In4);
[~,col] = size(x);
dum = zeros(1,col);
han(count) = plot(x,dum,'d-','DisplayName',['Catch trials,','N - ',num2str(sum([trials.present]==0))]);
hold on
In4 = ([trials.present]==0 & [trials.resp]==0);
x = find(In4);
[~,col] = size(x);
dum = zeros(1,col);
plot(x,dum,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color);


% set(gca,'box','off','tickdir','out')
% legend(han)
% xlabel('Trial No.')
% ylabel('contrast')
% title('Instance PEST iteration over trials')
% hold off

sess=[];
count = 1;
if (read_Data == 1)
    for i=1:numel(files)
        if i==numel(files)
            break
        end
        if ~(strcmp(files(i).name(1:13),files(i+1).name(1:13)))
            sess = [sess; i+1];
            figure(1);
            hold on
            line([sess(count) sess(count)],[0 max(contrast)],'Color',[0 1 1])
            count = count+1;
        end
        
%         if (files(i).name(end-6:end-4)=='001')
%             sess = [sess; i-1];  
%             figure(1);
%             hold on  
%             line([sess(count) sess(count)],[0 max(contrast)],'Color',[0 1 1])
%             count = count+1;
%         end
    end
else    
    for i=1:numel(filenames)
        if i==numel(filenames)
            break
        end
        if ~(strcmp(filenames(i).name(1:13),filenames(i+1).name(1:13)))
            sess = [sess; i+1];
            figure(1);
            hold on
            line([sess(count) sess(count)],[0 max(contrast)],'Color',[0 1 1])
            count = count+1;
        end
%         if (filenames(i).name(end-6:end-4)=='001')
%             sess = [sess; i-1];  
%             figure(1);
%             hold on  
%             line([sess(count) sess(count)],[0 max(contrast)],'Color',[0 1 1])
%             count = count+1;
%         end
    end
end
set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend(han)
xlabel('Trial No.')
% ylabel('contrast')
ylabel('Contrast')
title('Iteration of PEST levels')
hold off
saveas(figure(1),[sub_fig_path '\PEST.png']);
saveas(figure(1),[sub_fig_path '\PEST.fig']);
saveas(figure(1),[sub_fig_path '\PEST'],'epsc');