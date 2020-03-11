path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
% subject = 'A013';
% subject = 'Nikunj';
% subject = 'A092';
subject = 'A036';


pathtodata = fullfile(path,subject);

fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures_COMBINED';

sub_fig_path = fullfile(fig_path,subject);
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


path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\TIME_EXP_DATA';
pathtodata = fullfile(path,subject);

read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt1] = preprocessing(data); % this is a standard preprocessing call

[valid1,counter1] = countingTrialsNK(ppt1);

%% For Con. changing exp
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

%% For Time changing exp
trials1 = struct;
count = 1;
for i=1:length(ppt1)
%     if valid.validTrials(i)
    if valid1.drift(i)
        trials1(count).contrast = ppt1{i,1}.contrast;
%         trials(count).contrast = Untitled3(ppt1{i,1}.contrast,ppt1{i,1}.eccentricity,ppt1{i,1}.spatialFreq);
        trials1(count).resp = ppt1{i,1}.resp;
        trials1(count).present = ppt1{i,1}.present;
        trials1(count).respTime = ppt1{i,1}.responseTime;
        trials1(count).ecc = ppt1{i,1}.eccentricity;
        trials1(count).presTime = ppt1{i,1}.presTime;
        trials1(count).spFreq = ppt1{i,1}.spatialFreq;

        trials1(count).pixelAngle = ppt1{i,1}.pixelAngle;
        trials1(count).xoffset = ppt1{i,1}.xoffset;
        trials1(count).yoffset = ppt1{i,1}.yoffset;
        trials1(count).cuePrf = ppt1{i,1}.cuePrf;
        trials1(count).fixOn = ppt1{i,1}.fixOn;
        trials1(count).presTime = ppt1{i,1}.presTime;

        trials1(count).cueCtr = ppt1{i,1}.cueCtr;
        trials1(count).saccOn = ppt1{i,1}.saccOn;
        trials1(count).saccOff = ppt1{i,1}.saccOff;
        trials1(count).flashOn = ppt1{i,1}.flashOn;
        trials1(count).rampOff = ppt1{i,1}.rampOff;
        trials1(count).stimOff = ppt1{i,1}.stimOff;
        trials1(count).quit = ppt1{i,1}.quit; 
        
        %Correct for erroenous PEST levels
        if trials1(count).contrast == 0 && trials1(count).present == 1
            trials1(count).present = 0;
        end
%         if trials1(count).contrast ~= 0 && trials1(count).present == 0
%             trials1(count).contrast = 0;
%         end
        
        if trials1(count).contrast > 0.5
            trials1(count).contrast = 0.5;
        end
        count = count+1;
    end
end

ecc_levels1 = unique([trials1.ecc]);
contrast1 = [trials1.contrast];
pres_Time1 = [trials1.presTime];
sp_Freq1 = [trials1.spFreq];

count = 1;
for i=1:length(ecc_levels1)
    In1=([trials1.ecc]==ecc_levels1(i)); %indices%    
%     pres_levels= unique(pres_Time(In1));
    spFreq_lev = unique(sp_Freq1(In1));
    [~,t] = size(spFreq_lev);
%     for p = 1:length(pres_levels) %different presentation times   
    for j = 1:t %different spatial frequencies
%         In2=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i)); %indices
        In2=([trials1.spFreq]==spFreq_lev(j) & [trials1.ecc]==ecc_levels1(i)); %indices
%         freq_levels = unique(sp_Freq(In2));
%         for c=1:length(freq_levels)
            plotting1(count).ecc = ecc_levels1(i);
%             plotting(count).presTime = pres_levels(p);
            plotting1(count).spFreq = spFreq_lev(j);
%             In3=([trials.presTime]==pres_levels(p) & [trials.ecc]==ecc_levels(i) & [trials.spFreq]==freq_levels(c) & [trials.present]==1); %indices
            In3_present=([trials1.ecc]==ecc_levels1(i) & [trials1.spFreq]==spFreq_lev(j) & [trials1.present]==1); %indices
%             con_lev= unique(contrast(In3)); 
            pres_levels= unique(pres_Time1(In3_present)); 
            plotting1(count).presTime = pres_levels;
            x_present = find(In3_present); %Convert Logical into Sequential Index Vector
            plotting1(count).contrast = trials1(x_present(end)).contrast; %ALL TRIALS IN A CONDITION HAVE SAME CONTRAST            
%             for k = 1:length(con_lev)
            for k = 1:length(pres_levels)
                noTrials = sum([trials1.presTime]==pres_levels(k) & [trials1.ecc]==ecc_levels1(i)& [trials1.spFreq]==spFreq_lev(j));
                noStim = sum([trials1.present]==1 & [trials1.presTime]==pres_levels(k) & [trials1.ecc]==ecc_levels1(i) & [trials1.spFreq]==spFreq_lev(j));%stimulus displayed
                hits = sum([trials1.present]==1 &[trials1.resp]==1 & [trials1.presTime]==pres_levels(k) & [trials1.ecc]==ecc_levels1(i) & [trials1.spFreq]==spFreq_lev(j));
                misses = sum([trials1.present]==1 &[trials1.resp]==0 & [trials1.presTime]==pres_levels(k) & [trials1.ecc]==ecc_levels1(i) & [trials1.spFreq]==spFreq_lev(j));
                
                falseAlarms = sum([trials1.present]==0 &[trials1.resp]==1 & [trials1.presTime]==pres_levels(k));
                lures = sum([trials1.present]==0 & [trials1.presTime]==pres_levels(k));
                
                plotting1(count).noOfTrials(k) = noTrials; 
                plotting1(count).noOfStimDisp(k) = noStim;
                plotting1(count).hits(k) = hits;
                plotting1(count).misses(k) = misses;
                plotting1(count).percentCorrect(k) = 100 * (hits/noStim) ;
                
                plotting1(count).falseAlarms(k) = falseAlarms;
                plotting1(count).lures(k) = lures;                
                [plotting1(count).dprime(k), ~] = dprime(hits,falseAlarms,noStim,lures);                
            end  
            count = count+1;
%         end
    end
end

%% Presentation time vs Contrast
fh=1;
for i=1:length(ecc_levels1)
    for j=1:length(spFreq_lev)
%         In1=([plotting1.ecc]==ecc_levels1(i) & [plotting1.spFreq]==spFreq_lev(j)); %indices%    
%         p_times = plotting1(In1).presTime;
%         single_con = plotting1(In1).contrast;
%         
        In=([plotting.ecc]==ecc_levels1(i) & [plotting.spFreq]==spFreq_lev(j)); %indices% CON_EXP
       
        rows = find(In);
%         rows1 = find(In1);
%         range = max(plotting1(In1).presTime);%Max time
%         if range<500
%             range=500;
%         end
        
        for r1=1:length(rows)            
            p_times = plotting(rows(r1)).presTime;%ONLY Getting one values
            contrasts = plotting(rows(r1)).contrast;
            figure(fh)
%             scatter3(contrasts,repmat(p_times, 1,length(contrasts)),plotting(rows(r1)).percentCorrect,...
%                 'filled','DisplayName',['Fixed Pres Time - ',num2str(p_times)])
            plot3(contrasts,repmat(p_times, 1,length(contrasts)),plotting(rows(r1)).percentCorrect,...
                '-o','DisplayName',['Fixed Pres Time - ',num2str(p_times)])
            
%             imagesc(contrasts,repmat(p_times, 1,length(contrasts)),repmat(2,9))
%             imagesc(contrasts,repmat(p_times, 1,length(contrasts)),repmat(plotting(rows(r1)).percentCorrect,length(contrasts)))
%             pcolor(contrasts,repmat(p_times, 1,length(contrasts)),repmat(plotting(rows(r1)).percentCorrect,length(contrasts)))
            hold on
%             pause
        end
        
        hold on
        In1=([plotting1.ecc]==ecc_levels1(i) & [plotting1.spFreq]==spFreq_lev(j)); %indices% TIME_EXP
        %         rows1 = find(In1);
        p_times1 = plotting1(In1).presTime;%ONLY Getting one values
        single_contrast = plotting1(In1).contrast;
%         figure(fh)
%         scatter3(repmat(single_contrast,1,length(p_times1)),p_times1,plotting1(In1).percentCorrect,...
%             'filled','DisplayName',['Fixed Contrast - ',num2str(single_contrast)])
        plot3(repmat(single_contrast,1,length(p_times1)),p_times1,plotting1(In1).percentCorrect,...
            '-o','DisplayName',['Fixed Contrast - ',num2str(single_contrast)])
            
        xlabel('Contrast')
        ylabel('Presentation Time')
        zlabel('Proportion Reported Present')
        set(gca,'box','off','tickdir','out','FontSize',20)
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        legend
        title(['Ecc. - ',num2str(ecc_levels1(i)),',Freq. - ',num2str(spFreq_lev(j))])
%         p_times1 = unique(plotting(In).presTime);%ONLY Getting one values
%         single_con = plotting(In).contrast;
        
%         imagesc(single_con,p_times1,plotting(In).percentCorrect)
        
%         imagesc(single_con,repmat(p_times1, 1,length(single_con)),repmat(2,9))
        saveas(figure(fh),[sub_fig_path '/Ecc-',num2str(ecc_levels1(i)),'_Freq-',num2str(spFreq_lev(j))],'fig');
        saveas(figure(fh),[sub_fig_path '/Ecc-',num2str(ecc_levels1(i)),'_Freq-',num2str(spFreq_lev(j))],'epsc');
        fh=fh+1;
    end
end
