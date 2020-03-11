% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj';

path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
% subject = 'A013';%'Nikunj'
% subject = 'Nikunj';
% subject = 'A092';
subject = 'A036';
pathtodata = fullfile(path,subject);

read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
speed_var = 1;
if (read_Data ==1)
    if speed_var==1
        [data, files] = readdata(pathtodata, CalList_Vel_2(), true); %set to true to save in mat file and then use load. 
    else
        [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load.
    end
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

[valid,counter] = countingTrialsNK(ppt);

%if using Matlab 2017a only then will the following function work properly
% displayTrial(ppt);
% return

[noFiles,~] = size(ppt);
trials = struct;
% displayTrialSequence_nk(ppt,noFiles);
%      displayTrialSequence_nk(ppt,*array of numbers*);
%%
misTagTr = find(valid.misTags);
for i=1:length(misTagTr)
    inp = misTagTr(i);
% for i=1:length(ppt)
%     prompt = ['Enter a trial no between 1 and ',num2str(length(ppt)),' ,Input - '];
%     inp = input(prompt);
    
    if isempty(inp)
        break
%     elseif strcmp(num2str(inp),'d')
%         inp = inp+1;    
    end
%     inp = i;
    Quit = ppt{inp}.quit;
    CueCtr = ppt{inp}.cueCtr;
    StimOFF = ppt{inp}.stimOff;
    SaccON = ppt{inp}.saccOn;
    SaccOFF = ppt{inp}.saccOff;
    figure(1)
%     subplot(2,1,1);
    hold on
    X = double(ppt{inp}.x.position) + ppt{inp}.xoffset * ppt{inp}.pixelAngle;
    Y = double(ppt{inp}.y.position) + ppt{inp}.yoffset * ppt{inp}.pixelAngle;
    
    hx =  plot(X, 'Color', [0 0 200] / 255, 'HitTest', 'off'); % blue
    hy =  plot(Y, 'Color', [0 180 0] / 255, 'HitTest', 'off'); % green
    
    [~,emat] = (min(abs(ppt{inp}.saccades.start - SaccON)));
    if isempty(emat)
        continue
    end
    EMAT_saccOn = ppt{inp}.saccades.start(emat);
    EMAT_saccOff = EMAT_saccOn + ppt{inp}.saccades.duration(emat);
    h_EYERIS_saccOn = line([SaccON SaccON],[-200 500],'Color','green','LineWidth', 2);
    h_EMAT_saccOn = line([EMAT_saccOn EMAT_saccOn],[-200 500],'LineStyle','--','Color','green','LineWidth', 2);
    h_EMAT_saccOff = line([EMAT_saccOff EMAT_saccOff],[-200 500],'LineStyle','--','Color','blue','LineWidth', 2);    
    
    [~,emat] = (min(abs(ppt{inp}.drifts.start - SaccOFF)));
    if isempty(emat)
        continue
    end
    EMAT_driftOn = ppt{inp}.drifts.start(emat);
    EMAT_driftOff = EMAT_driftOn + ppt{inp}.drifts.duration(emat);
    h_EMAT_driftOn = line([EMAT_driftOn EMAT_driftOn],[-200 500],'LineStyle','--','Color','magenta','LineWidth', 2);
    h_EMAT_driftOff = line([EMAT_driftOff EMAT_driftOff],[-200 500],'LineStyle','--','Color','black','LineWidth', 2);    
    
%     h_stimOff = line([StimOFF StimOFF],[-200 500],'Color', [175 80 20]/255,'LineWidth', 2);
    h_stimOff = line([StimOFF StimOFF],[-200 500],'Color', 'red','LineWidth', 2);
    
    %Not displaying EYERIS saccOff as its close to the EMAT saccOff
%     h_EYERIS_saccOff = line([SaccOFF SaccOFF],[-200 500],'Color',[255 175 100]/255,'LineWidth', 2);
    h_EYERIS_saccOff = line([SaccOFF SaccOFF],[-200 500],'Color','blue','LineWidth', 2);
    leg = legend([hx hy h_EMAT_saccOn h_EYERIS_saccOn h_EYERIS_saccOff h_EMAT_saccOff h_EMAT_driftOn h_EMAT_driftOff h_stimOff],...
        'X', 'Y', 'EMAT_saccOn','EYERIS_saccOn','EYERIS_saccOff','EMAT_saccOff', 'EMAT_driftOn','EMAT_driftOff', 'StimOff',...
       'Location','EastOutside','Orientation','Vertical');
    
%     leg = legend([hx hy h_EMAT_saccOn h_EYERIS_saccOn h_EMAT_saccOff h_EMAT_driftOn h_EMAT_driftOff h_stimOff],...
%         'X', 'Y', 'EMAT_saccOn','EYERIS_saccOn','EMAT_saccOff', 'EMAT_driftOn','EMAT_driftOff', 'StimOff',...
%        'Location','EastOutside','Orientation','Vertical');
    set(leg,'Interpreter', 'none')
    legend boxoff
    grid on
    xlabel('time [ms]')
    ylabel('arcmin')
%     displayTrialSequence_nk_2(ppt,i);
%     displayTrialSequence_nk_2(ppt,inp);
%     legend('off')
    set(gca,'box','off','tickdir','out','FontSize',20)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    
    title(['Trial No. ',num2str(inp), 'Type - ', valid.name(inp)]);
%     subplot(2,1,2);
% %     tm = ((1:length(ppt{inp}.velocity)) .* (1/331) .* 1000);
%     plot(ppt{inp}.velocity);
%     hold on
%     if speed_var ==1
%          speedline_on = line([0 Quit],[ppt{inp}.velocity_on ppt{inp}.velocity_on],'Color',[255 175 100]/255,'LineWidth', 2,...
%              'DisplayName',['Saccade Speed On- ',num2str(ppt{inp}.velocity_on)]);
%          speedline_off = line([0 Quit],[ppt{inp}.velocity_off ppt{inp}.velocity_off],'Color',[0 175 100]/255,'LineWidth', 2,...
%              'DisplayName',['Saccade Speed Off- ',num2str(ppt{inp}.velocity_off)]);
%     else
%         speedline_on = line([0 Quit],[900 900],'Color',[255 175 100]/255,'LineWidth', 2,...
%              'DisplayName',['Saccade Speed On- ',num2str(900)]);
%         speedline_off = line([0 Quit],[300 300],'Color',[0 175 100]/255,'LineWidth', 2,...
%             'DisplayName',['Saccade Speed Off- ',num2str(300)]);
%     end
%     legend([speedline_on speedline_off],'Location','EastOutside','Orientation','Vertical');
%     legend boxoff
%     xlabel('time [ms]')
%     ylabel('Speed(arcmin/sec)')    
%     grid on
% %     xlim([CueCtr Quit]); 
% %     xlim([CueCtr (StimOFF+100)]);
%     title(['Trial No. ',num2str(inp)]);
%     set(gca,'box','off','tickdir','out','YScale', 'log')
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
%     waitforbuttonpress

    pause
    close all;
end
%%
% for i=1:length(ppt)
% %     trials.contrast(i) = ppt{i,1}.contrast;
% %     trials.resp(i) = ppt{i,1}.resp;
% %     trials.present(i) = ppt{i,1}.present;
%         trials(i).contrast = ppt{i,1}.contrast;
%         trials(i).resp = ppt{i,1}.resp;
%         trials(i).present = ppt{i,1}.present;
%         trials(i).respTime = ppt{i,1}.responseTime;
%         trials(i).ecc = ppt{i,1}.eccentricity;
%         trials(i).presTime = ppt{i,1}.presTime;
%         trials(i).spFreq = ppt{i,1}.spatialFreq;
% 
%         trials(i).pixelAngle = ppt{i,1}.pixelAngle;
%         trials(i).xoffset = ppt{i,1}.xoffset;
%         trials(i).yoffset = ppt{i,1}.yoffset;
%         trials(i).cuePrf = ppt{i,1}.cuePrf;
%         trials(i).fixOn = ppt{i,1}.fixOn;
%         trials(i).presTime = ppt{i,1}.presTime;
% 
%         trials(i).cueCtr = ppt{i,1}.cueCtr;
%         trials(i).saccOn = ppt{i,1}.saccOn;
%         trials(i).saccOff = ppt{i,1}.saccOff;
%         trials(i).flashOn = ppt{i,1}.flashOn;
%         trials(i).rampOff = ppt{i,1}.rampOff;
%         trials(i).stimOff = ppt{i,1}.stimOff;
%         trials(i).quit = ppt{i,1}.quit;         
%         
%         trials(i).saccade1 = vector2events(data.triggers{i}.saccade1);
%         trials(i).msaccade1 = vector2events(data.triggers{i}.msaccade1);
%         trials(i).em_event1 = vector2events(data.triggers{i}.em_event1);
%         
%         for( j = 1 : size(trials(i).em_event1.start,2) )        
%             for( k = 1 : size(trials(i).msaccade1.start,2) )
%                 if( trials(i).msaccade1.start(k) == trials(i).em_event1.start(j) + trials(i).em_event1.duration(j) )
%                     trials(i).msaccade1.start(k) = trials(i).em_event1.start(j);
%                     trials(i).msaccade1.duration(k) = trials(i).em_event1.duration(j);
%                 end
%             end
%             for( k = 1 : size(trials(i).saccade1.start,2) )
%                 if( trials(i).saccade1.start(k) == trials(i).em_event1.start(j) + trials(i).em_event1.duration(j) )
%                     trials(i).saccade1.start(k) = trials(i).em_event1.start(j);
%                     trials(i).saccade1.duration(k) = trials(i).em_event1.duration(j);
%                 end
%             end
%         end
% end


% for i=1:length(ppt)
% % for i=11
%     prompt = ['Enter a trial no between 1 and ',num2str(length(ppt)),' ,Input - '];
% %     x = input(prompt);
% %     displayTrialSequence_nk(ppt,x);
%     displayTrialSequence_nk(ppt,i);
% %      displayTrialSequence_nk(ppt,*array of numbers*);
%     title(['Trial No. ',num2str(x)]);
%     set(gca,'box','off','tickdir','out')
%     % Enlarge figure to full screen.
%     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% %     pause
% %     close all;
% end`