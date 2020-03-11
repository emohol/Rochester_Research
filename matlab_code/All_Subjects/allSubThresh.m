path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
% subject = 'A013';
% subject = 'Nikunj';
% subject = 'A092';
% subject = 'A036';
% pathtodata = fullfile(path,subject);
subjects = ["A013","Nikunj","A092","A036"];
% subjects = ["A013","Nikunj","A092"];
% fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\ALL_SUB';
fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\ALL_SUB\perf50';
% sub_fig_path = fullfile(fig_path,subject,'psychfit');
% sub_fig_path = fullfile(fig_path,subject,'land30');
if (~isfolder(fig_path))
    mkdir(fig_path)
end
avg_sens_2 = struct;
avg_sens_10 = struct;
load_mat =1;%set to 1 if loading from saved matfile
th_per = 0.5;%Threshold calculated at this performance
for si = 1:length(subjects)    
    if load_mat ==1
        matfilename = strcat(subjects(si),'_plot.mat');
        load(matfilename,'plotting')
        ecc_levels = unique([plotting.ecc]);
        contrast = [plotting.contrast];
        pres_Time = [plotting.presTime];
        sp_Freq = [plotting.spFreq];
    else        
        pathtodata = fullfile(path,subjects(si));
        bin = 1; %set to 1 if binning
        read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
        if (read_Data ==1)
            [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load.
        else
            load(fullfile(pathtodata, 'results.mat'));
        end
        [ppt] = preprocessing(data); % this is a standard preprocessing call
        
        [valid,counter] = countingTrialsNK(ppt);
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
        plotting = struct;
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
                    count = count+1;
                end
            end
        end
    end
    
    for i=1:12
        
        lvl = plotting(i).contrast; % contrast
        hits = plotting(i).hits; % # correct
        %     tr = plotting(i).noOfTrials;
        tr = plotting(i).noOfStimDisp;
        fa_rate = plotting(i).falseAlarms(1)/plotting(i).lures(1);
        
        [thresh, par]=psyfit(lvl,hits,tr,'Title', ['Ecc. - ',num2str(plotting(i).ecc)...
        ,' ,Pres. - ',num2str(plotting(i).presTime),' ,Freq. - ',num2str(plotting(i).spFreq)],'PlotOff','Extra',...
        'Chance',fa_rate,'Lapses','Auto','Thresh',th_per);
        plotting(i).thresh = round(thresh,2);%Round because that's the resolution of the shader       
        
        
    end
    
%     avg_sens_10 = [];
    figure(2)
    hold on
    for i=1:length(ecc_levels)
        In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
        In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
        t1 = plotting(In1).thresh;
        t2 = plotting(In2).thresh;
        std1 = plotting(In1).std;
        std2 = plotting(In2).std;
        stddev = [std1 std2];
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
        plotting(In1).sens = 1/t1;
        plotting(In2).sens = 1/t2;
        avg_sens_10(si).ecc(i) = ecc_levels(i);
        avg_sens_10(si).pres1(i) = 50;
        avg_sens_10(si).sens1(i) = thresh(1);
        avg_sens_10(si).pres2(i) = 500;
        avg_sens_10(si).sens2(i) = thresh(2);
%         avg_sens_10 = [avg_sens_10;sens];
        
        %             if sens>0.5
        %                 sens = 0.5;
        %             elseif sens<0
        %                 sens = 0;
        %             end
        if si == 1
            
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'s','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'s','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'s','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'s','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'s','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'s','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end
            
        elseif si == 2
            
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'v','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'v','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'v','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'v','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'v','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'v','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end
        elseif si==3
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'o','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'o','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'o','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'o','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'o','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'o','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end
         elseif si==4
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'x','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'x','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - A036,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'x','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'x','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - A036,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'x','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'x','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - A036,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end    
            
        end
    end
    
    hold off


    figure(1)
    hold on
    for i=1:length(ecc_levels)
        In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
        In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
        t1 = plotting(In1).thresh;
        t2 = plotting(In2).thresh;
        std1 = plotting(In1).std;
        std2 = plotting(In2).std; 
%         if std1==0
        stddev = [std1 std2];
        
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
        plotting(In1).sens = 1/t1;
        plotting(In2).sens = 1/t2;
        avg_sens_2(si).ecc(i) = ecc_levels(i);
        avg_sens_2(si).pres1(i) = 50;
        avg_sens_2(si).sens1(i) = thresh(1);
        avg_sens_2(si).pres2(i) = 500;
        avg_sens_2(si).sens2(i) = thresh(2);
        %             if sens>0.5
        %                 sens = 0.5;
        %             elseif sens<0
        %                 sens = 0;
        %             end
        if si == 1
            
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'s','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'s','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'s','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'s','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'s','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'s','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - A013,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end
            
        elseif si == 2
            
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'v','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'v','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'v','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'v','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'v','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'v','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - Nikunj,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end
        elseif si==3
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'o','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'o','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'o','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'o','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'o','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'o','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end
         elseif si==4
            if ecc_levels(i)==0
%                 errorbar([50 500],sens,stddev,'x','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'x','linewidth', 2,'Color','k','MarkerEdgeColor','k','MarkerSize', 15,...
                    'DisplayName',['Sub - A036,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==4
%                 errorbar([50 500],sens,stddev,'x','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'x','linewidth', 2,'Color','r','MarkerEdgeColor','r','MarkerSize', 15,...
                    'DisplayName',['Sub - A036,','Ecc - ',num2str(ecc_levels(i))])
                hold on
            elseif ecc_levels(i)==8
%                 errorbar([50 500],sens,stddev,'x','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
%                     'DisplayName',['Sub - A092,','Ecc - ',num2str(ecc_levels(i))])
                plot([50 500],thresh,'x','linewidth', 2,'Color','b','MarkerEdgeColor','b','MarkerSize', 15,...
                    'DisplayName',['Sub - A036,','Ecc - ',num2str(ecc_levels(i))])
                hold on

            end    
            
        end
    end   
    
    hold off   
    
end

figure(2)
hold on
avg_0_50 = mean([avg_sens_10(1).sens1(1) avg_sens_10(2).sens1(1) avg_sens_10(3).sens1(1) avg_sens_10(4).sens1(1)]);%Ecc 0, Pres 50
avg_4_50 = mean([avg_sens_10(1).sens1(2) avg_sens_10(2).sens1(2) avg_sens_10(3).sens1(2)  avg_sens_10(4).sens1(2)]);%Ecc 0, Pres 50
avg_8_50 = mean([avg_sens_10(1).sens1(3) avg_sens_10(2).sens1(3) avg_sens_10(3).sens1(3) avg_sens_10(4).sens1(3)]);%Ecc 0, Pres 50
avg_0_500 = mean([avg_sens_10(1).sens2(1) avg_sens_10(2).sens2(1) avg_sens_10(3).sens2(1) avg_sens_10(4).sens2(1)]);%Ecc 0, Pres 500
avg_4_500 = mean([avg_sens_10(1).sens2(2) avg_sens_10(2).sens2(2) avg_sens_10(3).sens2(2) avg_sens_10(4).sens2(2)]);%Ecc 0, Pres 500
avg_8_500 = mean([avg_sens_10(1).sens2(3) avg_sens_10(2).sens2(3) avg_sens_10(3).sens2(3) avg_sens_10(4).sens2(3)]);%Ecc 0, Pres 500
plot([50 500],[avg_0_50 avg_0_500],'-*','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize', 10,...
                    'DisplayName',['ALL - SUB,','Ecc - 0'])
hold on
plot([50 500],[avg_4_50 avg_4_500],'-*','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize', 10,...
                    'DisplayName',['ALL - SUB,','Ecc - 4'])
hold on
plot([50 500],[avg_8_50 avg_8_500],'-*','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize', 10,...
                    'DisplayName',['ALL - SUB,','Ecc - 8'])
set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('Location','westoutside')
xticks([50 500])
xlim([-50 600])
xlabel('Presentation time')
ylabel('Threshold')
title('High SpFreq')
saveas(figure(2),[fig_path '\Thresh_Pres_10'],'epsc');
hold off


figure(1)
hold on
avg_0_50 = mean([avg_sens_2(1).sens1(1) avg_sens_2(2).sens1(1) avg_sens_2(3).sens1(1) avg_sens_2(4).sens1(1)]);%Ecc 0, Pres 50
avg_4_50 = mean([avg_sens_2(1).sens1(2) avg_sens_2(2).sens1(2) avg_sens_2(3).sens1(2)  avg_sens_2(4).sens1(2)]);%Ecc 0, Pres 50
avg_8_50 = mean([avg_sens_2(1).sens1(3) avg_sens_2(2).sens1(3) avg_sens_2(3).sens1(3) avg_sens_2(4).sens1(3)]);%Ecc 0, Pres 50
avg_0_500 = mean([avg_sens_2(1).sens2(1) avg_sens_2(2).sens2(1) avg_sens_2(3).sens2(1) avg_sens_2(4).sens2(1)]);%Ecc 0, Pres 500
avg_4_500 = mean([avg_sens_2(1).sens2(2) avg_sens_2(2).sens2(2) avg_sens_2(3).sens2(2) avg_sens_2(4).sens2(2)]);%Ecc 0, Pres 500
avg_8_500 = mean([avg_sens_2(1).sens2(3) avg_sens_2(2).sens2(3) avg_sens_2(3).sens2(3) avg_sens_2(4).sens2(3)]);%Ecc 0, Pres 500
plot([50 500],[avg_0_50 avg_0_500],'-*','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize', 10,...
                    'DisplayName',['ALL - SUB,','Ecc - 0'])
hold on
plot([50 500],[avg_4_50 avg_4_500],'-*','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize', 10,...
                    'DisplayName',['ALL - SUB,','Ecc - 4'])
hold on
plot([50 500],[avg_8_50 avg_8_500],'-*','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize', 10,...
                    'DisplayName',['ALL - SUB,','Ecc - 8'])
set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('Location','westoutside')
xticks([50 500])
xlim([-50 600])
xlabel('Presentation time')
ylabel('Threshold')
title('Low SpFreq')
saveas(figure(1),[fig_path '\Thresh_Pres_2'],'epsc');