function [valid,counter] = countingTrialsNK_DDPI(ppt,sRate)%,numTooShort,numBlinks,numNoTracks,numSaccades,numMicroSaccades)
% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj_NoMask';
% % pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj';
% read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
% if (read_Data ==1)
%     [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
% else
%     load(fullfile(pathtodata, 'results.mat'));
% end
% [ppt] = preprocessing(data); % this is a standard preprocessing call
numTotalTrials = length(ppt);

% ppt = osRemoverInEM(ppt);

% counter.TooShort = 0;
% counter.BigOffset = 0;
counter.Blinks = 0;
counter.NoTracks = 0;
% counter.NoResponse = 0;
% counter.manualDiscard = 0;
counter.Saccades = 0;
% counter.Fixation = 0;
% counter.Span = 0;
counter.Microsaccades = 0;
counter.Drifts = 0;
counter.manualDiscard = 0;
counter.validTrials =0;
counter.noTags =0;
counter.totalTrials = numTotalTrials;
valid = struct(...;
    'name',strings(size(ppt)),...;
    'peakVelocity',zeros(size(ppt),'double'),...;
    'IndexPeakVelocity',zeros(size(ppt)),...;
    'noTags', false(size(ppt)),...;
    'blink', false(size(ppt)),...;
    'notrack', false(size(ppt)),...;
    'sac', false(size(ppt)),...;
    'micsac', false(size(ppt)),...;
    'manualdiscard', false(size(ppt)),...;
    'validTrials', false(size(ppt)),...;
    'drift', false(size(ppt)));


% intersectEvents(Events, Start, Duration)

for ii = 1:length(ppt)
%     if ppt{ii}.present ==0
%         continue;
%     end
%     FixationON = ppt{ii}.fixOn;
%     CueCtr = ppt{ii}.cueCtr;  
    SaccON = ppt{ii}.saccOn * (sRate/1000);
%     FlashON = ppt{ii}.flashOn;
    SaccOFF = ppt{ii}.saccOff* (sRate/1000); 
    StimOFF = ppt{ii}.stimOff* (sRate/1000);
%     Quit = ppt{ii}.quit;
%     RespTime = ppt{ii}.responseTime;
%     timeOn = round(FlashON);

%     timeOn = round(CueCtr);%When the fixation dot disappears

    timeOn = round(SaccON);%Cause presentation time starts after central fixation
    timeOff = round(SaccOFF);
    
%     timeOn = round(FlashON);
%     timeOff = round(StimOFF);
%     valid.peakVelocity(ii)=max(ppt{ii}.velocity(round(CueCtr:end))); 
%     valid.IndexPeakVelocity(ii)=find(ppt{ii}.velocity == valid.peakVelocity(ii));
     
    if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.blinks)%%%BLINKS
        valid.blink(ii) = true;
        counter.Blinks = counter.Blinks+1;
        valid.name(ii) = 'Blink';
        continue;
    end
    
    if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.notracks)%%%NO TRACK TRIALS
        valid.notrack(ii) = true;
        counter.NoTracks = counter.NoTracks+1;
        valid.name(ii) = 'NoTrack';
        continue;
    end
    
    if  ~isempty(ppt{ii}.invalid.start)
        valid.manualdiscard(ii) = true;
        counter.manualDiscard = counter.manualDiscard+1;
        valid.name(ii) = 'ManualDiscard';
        continue;
    end
    
%     if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS
    if isIncludedInPeriod(ppt{ii}.drifts,timeOn,timeOff-timeOn)%%%DRIFT TRIALS
        %     if isIntersectedIn(timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS
        valid.drift(ii) = true ;
        counter.Drifts = counter.Drifts+1;
        valid.name(ii) = 'Drift';
        continue;
    end
    
    if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.saccades)%%%SACCADE TRIALS
        valid.sac(ii) = true;
        counter.Saccades = counter.Saccades+1;
        timeOn = round(SaccOFF);
        timeOff = round(StimOFF);
        valid.name(ii) = 'Saccade';
             if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS
%              if (isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.drifts) && ...
%                    isIncludedInPeriod(ppt{ii}.saccades,timeOn,timeOff-timeOn))  
        %     if isIntersectedIn(timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS
                valid.validTrials(ii) = true ;
                counter.validTrials = counter.validTrials+1;
                counter.Saccades = counter.Saccades-1;
                valid.name(ii) = 'VALID';
                continue;
            end
        continue;
    end  
    
    if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.microsaccades)%%%MICROSACCADE TRIALS
        valid.micsac(ii) = true;
        counter.Microsaccades = counter.Microsaccades+1;
        valid.name(ii) = 'MicroSaccade';
        continue;
    end
    
    valid.noTags(ii) = true;
    counter.noTags = counter.noTags+1;   
   
end

fprintf('\nThere are %i in total trials.\n', numTotalTrials)
fprintf('There are %i trials that have blinks.\n',counter.Blinks)
fprintf('There are %i trials do not track.\n',counter.NoTracks)
fprintf('There are %i trials that have been manually discarded.\n',counter.manualDiscard)
fprintf('There are %i trials that have only saccades in the stimulus presentation period.\n',counter.Saccades)
fprintf('There are %i trials that have only microsaccades in the stimulus presentation period.\n',counter.Microsaccades)
fprintf('There are %i trials that have drift in the saccade period.\n',counter.Drifts)
fprintf('There are %i trials that have no tagged events\n',counter.noTags)
fprintf('There are %i trials that are valid.\n',counter.validTrials)
% fprintf('There are %i trials after filter.\n',  sum(cell2mat(struct2cell(counter))) - counter.TotalTask)