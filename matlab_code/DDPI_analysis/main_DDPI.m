% pathtodata = pwd;
% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj_NoMask';
% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj';
pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\DDPI\MAC';
% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\test1';
% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\photodiode';
bin = 0; %set to 1 if binning
read_Data = 1; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
% [ppt] = preprocessing_DDPI(data); % this is a standard preprocessing call
[ppt] = preprocessing(data); % this is a standard preprocessing call
sRate = 331;
[valid,counter] = countingTrialsNK_DDPI(ppt,sRate);

%%
% %To display table in figure
% f = figure;
% uit = uitable(f);
% d={'Total Trials',counter.totalTrials; 'Blinks', counter.Blinks;...
%     'No Tracks',counter.NoTracks; 'Saccades only', counter.Saccades;...
%     'MicroSaccades only', counter.Microsaccades; 'Drifts only',counter.Drifts;...
%     'Discarded',(counter.manualDiscard+counter.weird); 'Valid Trials',counter.validTrials};
% uit.Data = d;
% uit.ColumnName = {};   
% uit.RowName = {};
% uit.Position = [150 150 250 250];
% uit.FontSize = 10; 

% %To save figure to csv to load in Latex Report
% rownames = {'Total';'Blink';'No Track';'Saccade only';'MicroSaccade only';'Drift only';'Discarded';'Valid'};
% col_1 = [counter.totalTrials;counter.Blinks;counter.NoTracks;counter.Saccades;counter.Microsaccades;counter.Drifts;...
%     (counter.manualDiscard+counter.weird);counter.validTrials];
% val = table(col_1,'RowNames',rownames,'VariableNames',{'Trials'});
% disp(val);
% writetable(val,'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\trials.csv','WriteRowNames',true,'WriteVariableNames',true);

trials = struct;
count = 1;
for i=1:length(ppt)
%     trials.contrast(i) = ppt{i,1}.contrast;
%     trials.resp(i) = ppt{i,1}.resp;
%     trials.present(i) = ppt{i,1}.present;
    if valid.validTrials(i)
%         trials(count).contrast = ppt{i,1}.contrast;
        trials(i).contrast = (ppt{i,1}.contrast)*2;
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
        count = count+1;
    end
    %use the date and time info to create 'session', meaning in one
    %sitting. Also multiple sessions on the same day!
end

%round off
if (bin==1)   
    for i=1:length(trials)
        trials(i).contrast = round([trials(i).contrast],3);
    end
else
    bin = 0;
end

% boxplot(x,{H,M})
% grpstats(x,{H,M})
% [u g]=grpstats(x,{H,M},{'mean','gname'})
% [u g n]=grpstats(x,{H,M},{'mean','gname','numel'})

ecc_levels = unique([trials.ecc]);
contrast = [trials.contrast];
pres_Time = [trials.presTime];
sp_Freq = [trials.spFreq];

% x=sort([trials.contrast]);
T = table;
count = 1;
for i=1:length(ecc_levels)
%     plotting = struct([]);
% for i=length(ecc_levels)
% for i=1:1:2
% for i=2:1:3  
%     count1=1;
    In1=([trials.ecc]==ecc_levels(i)); %indices
%     con_lev= unique(contrast(In1));    
%     In2=([trials.presTime]==pres_levels(i)); %indices
    pres_levels= unique(pres_Time(In1));
%     [~,u] = size(pres_levels);
    for p = 1:length(pres_levels) %different presentation times   
%         st = p; %to run the cor_rej plot only once per pres_level
        In2=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)); %indices
        freq_levels = unique(sp_Freq(In2));
        lures = sum([trials.present]==0 & [trials.presTime]==pres_levels(p));              
        cor_rej = sum([trials.present]==0 &[trials.resp]==0 & [trials.presTime]==pres_levels(p));
        falseAlarms = sum([trials.present]==0 &[trials.resp]==1 & [trials.presTime]==pres_levels(p));
%         figure(2*i); 
%         plot(0,100*(cor_rej/lures),'d','DisplayName',['0-Con.,Pres.-',num2str(pres_levels(p))]);
%         formatSpec = "(%d,%d)";
%         txt = sprintf(formatSpec,cor_rej,falseAlarms);
%         text(0,100*(cor_rej/lures),txt);
%         hold on
        for c=1:length(freq_levels)
            plotting(count).ecc = ecc_levels(i);
            plotting(count).presTime = pres_levels(p);
            plotting(count).spFreq = freq_levels(c);
            In3=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c)); %indices
            con_lev= unique(contrast(In3)); 
            plotting(count).contrast = con_lev;            
            for j = 1:length(con_lev)
                noTrials = sum([trials.contrast]==plotting(count).contrast(j) & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)& [trials.spFreq]==freq_levels(c));
                noStim = sum([trials.contrast] == plotting(count).contrast(j) & [trials.present]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));%stimulus displayed
                hits = sum([trials.contrast]==plotting(count).contrast(j) & [trials.present]==1 &[trials.resp]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));
                misses = sum([trials.contrast]==plotting(count).contrast(j) & [trials.present]==1 &[trials.resp]==0 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));
%                 falseAlarms = sum([trials.present]==0 &[trials.resp]==1 & [trials.presTime]==pres_levels(p)); % both lures and falseAlarms need to be based on your
                %'session', later on you want to do the same for hits and
                %no.StimulusDisplayed as well.
%                  falseAlarms = sum([trials.present]==0 &[trials.resp]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));
%                 lures = sum([trials.present]==0 & [trials.presTime]==pres_levels(p));
                plotting(count).noOfTrials(j) = noTrials; 
                plotting(count).noOfStimDisp(j) = noStim;
                plotting(count).hits(j) = hits;
                plotting(count).misses(j) = misses;
                plotting(count).percentCorrect(j) = 100 * (hits/noStim) ;
                plotting(count).dprime(j) = dprime(hits,falseAlarms,noStim,lures);
            end
%             figure(2*i - 1);
%             plot(plotting(count).contrast,plotting(count).dprime,'d-','DisplayName',['Pres.-',num2str(pres_levels(p)),',Freq.-',num2str(freq_levels(c))])
%     %         plot(plotting(count).contrast,plotting(count).dprime,'d-','DisplayName','Pres. time %d',pres_levels(p)) 
%             for j = 1:length(con_lev)
%                 formatSpec = "(%d,%d)";
%                 txt = sprintf(formatSpec,plotting(count).hits(j),plotting(count).misses(j));
%                 text(plotting(count).contrast(j),plotting(count).dprime(j),txt);
%             end
% %         if i==length(ecc_levels)
% %             if p==length(pres_levels)
%             if (c==length(freq_levels) && p==length(pres_levels))
%                 figure(2*i - 1);
%                 set(gca,'box','off','tickdir','out')
%                 % Enlarge figure to full screen.
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     %             for l=1:length(pres_levels)
%     %                 str = sprintf('Pres. time %d',pres_levels(l));
%     %                 legend(str)
%     %             end
%                 % legend('ecc 4','ecc 8')
%                 legend
%                 xlabel('contrast')
%                 ylabel('d-prime')
%                 title(['Eccentricity ' num2str(ecc_levels(i))])
%                 hold off
%                 if (bin==1)
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/binned/ecc_%d_d_prime.png';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i-1),txt1);
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/binned/ecc_%d_d_prime.fig';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i-1),txt1);
%                 else
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/unbinned/ecc_%d_d_prime.png';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i-1),txt1);
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/ecc_%d_d_prime.fig';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i-1),txt1);
%                 end
%             else
%                 hold on
%             end
        
            %Percent plots
%             figure(2*i);            
%             plot(plotting(count).contrast,plotting(count).percentCorrect,'d-','DisplayName',['Pres.-',num2str(pres_levels(p)),',Freq.-',num2str(freq_levels(c))]);
%     %         plot(plotting(count).contrast,plotting(count).dprime,'d-','DisplayName','Pres. time %d',pres_levels(p)) 
%             for j = 1:length(con_lev)
%                 formatSpec = "(%d,%d)";
%                 txt = sprintf(formatSpec,plotting(count).hits(j),plotting(count).misses(j));
%                 text(plotting(count).contrast(j),plotting(count).percentCorrect(j),txt);
%             end
%             if (c==length(freq_levels) && p==length(pres_levels))
%                 figure(2*i);
%                 set(gca,'box','off','tickdir','out')
%                 % Enlarge figure to full screen.
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%                 legend
%                 xlabel('contrast')
%                 ylabel('percentage correct')
%                 title(['Eccentricity ' num2str(ecc_levels(i))])
%                 hold off
%                 if (bin==1)
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/binned/ecc_%d_correct.png';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i),txt1);
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/binned/ecc_%d_correct.fig';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i),txt1);
%                 else
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/unbinned/ecc_%d_correct.png';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i),txt1);
%                     formatSpec1 = 'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/unbinned/ecc_%d_correct.fig';
%                     txt1 = sprintf(formatSpec1,ecc_levels(i));
%                     saveas(figure(2*i),txt1);
%                 end
%             else
%                 hold on
%             end
                       
            %PEST plot        
            In3_resp1=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1 &[trials.resp]==1); %indices
            con_lev_resp1 = contrast(In3_resp1);
            x_resp1 = find(In3_resp1); %Convert Logical into Sequential Index Vector
            
            In3_resp0=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1); %indices
            con_lev_resp0 = contrast(In3_resp0);
            x_resp0 = find(In3_resp0); %Convert Logical into Sequential Index Vector 
            last_con=trials(x_resp0(end)).contrast;
            
            figure(2*length(ecc_levels)+2);
            if (freq_levels(c) == 2 && ecc_levels(i)== 0)
                han(count) = plot(x_resp0,con_lev_resp0,':d','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1))...
                    ,',C - ',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (freq_levels(c) == 2 && ecc_levels(i)== 4)
                han(count) = plot(x_resp0,con_lev_resp0,'--d','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (freq_levels(c) == 2 && ecc_levels(i)== 8)
                han(count) = plot(x_resp0,con_lev_resp0,'-d','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (freq_levels(c) == 10 && ecc_levels(i)== 0)
                han(count) = plot(x_resp0,con_lev_resp0,':s','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))...
                   ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1))...
                   ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'s','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color) 
            elseif (freq_levels(c) == 10 && ecc_levels(i)== 4)
                han(count) = plot(x_resp0,con_lev_resp0,'--s','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'s','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            elseif (freq_levels(c) == 10 && ecc_levels(i)== 8)
                han(count) = plot(x_resp0,con_lev_resp0,'-s','DisplayName',['Pres - ',num2str(pres_levels(p)),',Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))...
                    ,',N - ',num2str(sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1))...
                    ,',Con.-',num2str(round(last_con,3))]);
                hold on
                plot(x_resp1,con_lev_resp1,'s','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color)
            end
            
             %response time plots
            In4=([trials.present]==1 &[trials.resp]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c)); %indices
            resp = minus([trials.respTime],[trials.stimOff]);
            resp(:) = resp(:)/1000; %get in seconds
            plotting(count).respTime = mean(resp(In4));
%             plotting(count).stdDev = std(resp(In4));
            plotting(count).stdMean = std(resp(In4))/sqrt(length(resp(In4)));
    %         figure(2*i);
%             figure(2*length(ecc_levels)+1);
%     %         a = repmat(plotting(count).ecc,1,length(plotting(count).respTime));%to make x and y plots equal in length
%             e=errorbar(plotting(count).ecc,plotting(count).respTime,(plotting(count).stdMean),'DisplayName',['Ecc.-',num2str(ecc_levels(i)),',Pres.-',num2str(pres_levels(p)),',Freq.-',num2str(freq_levels(c))]);
%             e.Marker = '*';
%             e.Color = han(count).Color;
%             e.MarkerSize = 12;
%             if (i==length(ecc_levels) && c==length(freq_levels) && p==length(pres_levels))
%     %         if p==length(pres_levels)
%     %             figure(2*i);
%                 figure(2*length(ecc_levels)+1);
%                 set(gca,'box','off','tickdir','out')
%                 % Enlarge figure to full screen.
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     %             for l=1:length(pres_levels)
%     %                 str = sprintf('Pres. time %d',pres_levels(l));
%     %                 legend(str)
%     %             end
%                 % legend('ecc 4','ecc 8')
%                 legend
%                 xlabel('eccentricity')
%                 ylabel('Response time(sec)')
%                 title(['Ecc vs Response-time errorbar'])
%                 hold off                           
%                 saveas(figure(2*length(ecc_levels)+1),'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/errorbar.png');
%                 saveas(figure(2*length(ecc_levels)+1),'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/errorbar.fig');
%             else
%                 hold on
%             end 
            
            total_trials = sum([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c));
            total_hits = sum([trials.present]==1 &[trials.resp]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)& [trials.spFreq]==freq_levels(c));
            total_stim = sum([trials.present]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)& [trials.spFreq]==freq_levels(c));%stimulus displayed
            total_fa = sum([trials.present]==0 &[trials.resp]==1 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)& [trials.spFreq]==freq_levels(c));
            total_lures = sum([trials.present]==0 & [trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)& [trials.spFreq]==freq_levels(c));
            x_t = table(plotting(count).ecc, plotting(count).presTime,plotting(count).spFreq, total_trials, total_hits, total_stim, total_fa,total_lures);
            T = [T; x_t];   
            count = count+1;
        end
    end    
%     plotting = struct([]);
end
figure(2*length(ecc_levels)+1);
% e=errorbar(plotting(count).ecc,plotting(count).respTime,(plotting(count).stdMean),'DisplayName',['Ecc.-',num2str(ecc_levels(i)),',Pres.-',num2str(pres_levels(p)),',Freq.-',num2str(freq_levels(c))]);
var =1;
x=[plotting(var).ecc plotting(var+4).ecc plotting(var+8).ecc];
y=[plotting(var).respTime plotting(var+4).respTime plotting(var+8).respTime];
err=[plotting(var).stdMean plotting(var+4).stdMean plotting(var+8).stdMean];
errorbar(x,y,err,'-d','MarkerSize',12,'LineWidth',2,'DisplayName',['Pres.-',num2str(plotting(var).presTime),',Freq.-',num2str(plotting(var).spFreq)]);
hold on
var =2;
x=[plotting(var).ecc plotting(var+4).ecc plotting(var+8).ecc];
y=[plotting(var).respTime plotting(var+4).respTime plotting(var+8).respTime];
err=[plotting(var).stdMean plotting(var+4).stdMean plotting(var+8).stdMean];
errorbar(x,y,err,'-x','MarkerSize',12,'LineWidth',2,'DisplayName',['Pres.-',num2str(plotting(var).presTime),',Freq.-',num2str(plotting(var).spFreq)]);
hold on
var =3;
x=[plotting(var).ecc plotting(var+4).ecc plotting(var+8).ecc];
y=[plotting(var).respTime plotting(var+4).respTime plotting(var+8).respTime];
err=[plotting(var).stdMean plotting(var+4).stdMean plotting(var+8).stdMean];
errorbar(x,y,err,'--d','MarkerSize',12,'LineWidth',2,'DisplayName',['Pres.-',num2str(plotting(var).presTime),',Freq.-',num2str(plotting(var).spFreq)]);
hold on
var =4;
x=[plotting(var).ecc plotting(var+4).ecc plotting(var+8).ecc];
y=[plotting(var).respTime plotting(var+4).respTime plotting(var+8).respTime];
err=[plotting(var).stdMean plotting(var+4).stdMean plotting(var+8).stdMean];
errorbar(x,y,err,'--x','MarkerSize',12,'LineWidth',2,'DisplayName',['Pres.-',num2str(plotting(var).presTime),',Freq.-',num2str(plotting(var).spFreq)]);
% hold off
% e=errorbar(plotting(1:4:12).ecc,plotting(1:4:12).respTime,(plotting(1:4:12).stdMean),...
%     '-d','MarkerSize',12,'DisplayName',['Pres.-',num2str(plotting(1).presTime),',Freq.-',num2str(plotting(1).spFreq)]);
% % e.Color = han(count).Color;
% e=errorbar(plotting(2:4:12).ecc,plotting(2:4:12).respTime,(plotting(2:4:12).stdMean),...
%     '-x','MarkerSize',12,'DisplayName',['Pres.-',num2str(plotting(2).presTime),',Freq.-',num2str(plotting(2).spFreq)]);
% e=errorbar(plotting(3:4:12).ecc,plotting(3:4:12).respTime,(plotting(3:4:12).stdMean),...
%     '--d','MarkerSize',12,'DisplayName',['Pres.-',num2str(plotting(3).presTime),',Freq.-',num2str(plotting(3).spFreq)]);
% e=errorbar(plotting(4:4:12).ecc,plotting(4:4:12).respTime,(plotting(4:4:12).stdMean),...
%     '--x','MarkerSize',12,'DisplayName',['Pres.-',num2str(plotting(4).presTime),',Freq.-',num2str(plotting(4).spFreq)]);
xlim([-1, 9]);
set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend
xlabel('eccentricity')
ylabel('Response time(sec)')
title(['Ecc vs Response-time errorbar'])
hold off                   
saveas(figure(2*length(ecc_levels)+1),'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/errorbar.png');
saveas(figure(2*length(ecc_levels)+1),'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/errorbar.fig'); 
            
T.Properties.VariableNames = {'Ecc' 'PresTime' 'spFreq' 'Trials' 'Hits' 'NumberStimulus' 'FalseAlarms' 'Lures'};

In4 = ([trials.present]==0);
x = find(In4);
[~,col] = size(x);
dum = zeros(1,col);
figure(2*length(ecc_levels)+2);
han(count) = plot(x,dum,'d-','DisplayName',['Catch trials,','N - ',num2str(sum([trials.present]==0))]);
hold on
In4 = ([trials.present]==0 & [trials.resp]==0);
x = find(In4);
[~,col] = size(x);
dum = zeros(1,col);
plot(x,dum,'d','MarkerFaceColor',han(count).Color,'MarkerEdgeColor',han(count).Color);

%To draw vertical lines to mark sessions
% sess=[];
% count = 1;
% if (read_Data == 1)
%     for i=1:numel(files)
%         if (files(i).name(end-6:end-4)=='001')
%             sess = [sess; i-1];  
%             figure(2*length(ecc_levels)+2)
%             hold on  
%             line([sess(count) sess(count)],[0 max(contrast)],'Color',[0 1 1])
%             count = count+1;
%         end
%     end
% else
%     for i=1:numel(filenames)
%         if (filenames(i).name(end-6:end-4)=='001')
%             sess = [sess; i-1];  
%             figure(2*length(ecc_levels)+2)
%             hold on  
%             line([sess(count) sess(count)],[0 max(contrast)],'Color',[0 1 1])
%             count = count+1;
%         end
%     end
% end

set(gca,'box','off','tickdir','out')
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend(han)
xlabel('Trial No.')
ylabel('contrast')
title('Instance PEST iteration over trials')
hold off
saveas(figure(2*length(ecc_levels)+2),'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/PEST.png');
saveas(figure(2*length(ecc_levels)+2),'C:/Users/Ruccilab/Box/Vis_Dynamics/Figures/PEST.fig');


y = resp;
g1 = [trials(:).ecc];
g2 = [trials(:).presTime];
g3 = [trials(:).spFreq];
g4 = [trials(:).resp];
g5 = [trials(:).present];
[~,~,stats]=anovan(y,{g1,g2,g3,g4,g5},'model','interaction','varnames',{'Ecc','Pres','Freq','Response','Stimulus'});
% results = multcompare(stats,'Dimension',[1 2]);