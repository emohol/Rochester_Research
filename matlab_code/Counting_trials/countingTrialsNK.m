function [valid,counter] = countingTrialsNK(ppt)%,numTooShort,numBlinks,numNoTracks,numSaccades,numMicroSaccades)
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
counter.Saccades1 = 0;
counter.Saccades2 = 0;
counter.SaccadesMore = 0;
% counter.Fixation = 0;
% counter.Span = 0;
counter.Microsaccades = 0;
counter.Microsaccades1 = 0;
counter.Microsaccades2 = 0;
counter.MicrosaccadesMore = 0;
counter.Drifts = 0;
counter.Discard = 0;
counter.validTrials =0;
counter.noTags =0;
counter.misTags = 0;
counter.landDist = 0;
counter.totalTrials = numTotalTrials;
valid = struct(...;
    'name',strings(size(ppt)),...;
    'peakVelocity',zeros(size(ppt),'double'),...;
    'IndexPeakVelocity',zeros(size(ppt)),...;
    'noTags', false(size(ppt)),...;
    'misTags', false(size(ppt)),...;
    'landDist', false(size(ppt)),...;
    'blink', false(size(ppt)),...;
    'notrack', false(size(ppt)),...;
    'sac', false(size(ppt)),...;
    'sac1', false(size(ppt)),...;
    'sac2', false(size(ppt)),...;
    'sacMore', false(size(ppt)),...;
    'micsac', false(size(ppt)),...;
    'micsac1', false(size(ppt)),...;
    'micsac2', false(size(ppt)),...;
    'micsacMore', false(size(ppt)),...;
    'discard', false(size(ppt)),...;
    'validTrials', false(size(ppt)),...;
    'drift', false(size(ppt)));


% intersectEvents(Events, Start, Duration)

for ii = 1:length(ppt)
%     if ppt{ii}.present ==0
%         continue;
%     end
%     FixationON = ppt{ii}.fixOn;
%     CueCtr = ppt{ii}.cueCtr;  
    SaccON = ppt{ii}.saccOn;
%     FlashON = ppt{ii}.flashOn;
    SaccOFF = ppt{ii}.saccOff; 
    StimOFF = ppt{ii}.stimOff;
%     Quit = ppt{ii}.quit;
%     RespTime = ppt{ii}.responseTime;
%     timeOn = round(FlashON);

%     timeOn = round(CueCtr);%When the fixation dot disappears

    timeOn = round(SaccON);%Cause presentation time starts after central fixation
    timeOff = round(SaccOFF);
    
    X = double(ppt{ii}.x.position) + ppt{ii}.xoffset * ppt{ii}.pixelAngle;
    Y = double(ppt{ii}.y.position) + ppt{ii}.yoffset * ppt{ii}.pixelAngle;
%     timeOn = round(FlashON);
%     timeOff = round(StimOFF);
%     valid.peakVelocity(ii)=max(ppt{ii}.velocity(round(CueCtr:end))); 
%     valid.IndexPeakVelocity(ii)=find(ppt{ii}.velocity == valid.peakVelocity(ii));

    if ~isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.saccades)%%%NO tagged saccade event
        valid.noTags(ii) = true;
        counter.noTags = counter.noTags+1;
        valid.name(ii) = 'noTags';
        valid.discard(ii) = true;
        counter.Discard = counter.Discard+1;
        continue;
    end
    
    if isIntersectedIn(timeOn,StimOFF-timeOn,ppt{ii}.blinks)%%%BLINKS
        valid.blink(ii) = true;
        counter.Blinks = counter.Blinks+1;
        valid.name(ii) = 'Blink';
        valid.discard(ii) = true;
        counter.Discard = counter.Discard+1;
        continue;
    end
    
    if isIntersectedIn(timeOn,StimOFF-timeOn,ppt{ii}.notracks)%%%NO TRACK TRIALS
        valid.notrack(ii) = true;
        counter.NoTracks = counter.NoTracks+1;
        valid.name(ii) = 'NoTrack';
        valid.discard(ii) = true;
        counter.Discard = counter.Discard+1;
        continue;
    end
    
%     if  ~isempty(ppt{ii}.invalid.start)
%         valid.discard(ii) = true;
%         counter.Discard = counter.Discard+1;
%         valid.name(ii) = 'ManualDiscard';
%         continue;
%     end
    
    if isIncludedInPeriod(ppt{ii}.drifts,timeOn,timeOff-timeOn)
        valid.misTags(ii) = true;
        counter.misTags = counter.misTags+1;   
        valid.name(ii) = 'MisTags';
        valid.discard(ii) = true;
        counter.Discard = counter.Discard+1;
        continue;
    end
    
    %To calculate landing distance, at EYERIS-saccOff
    landDist = sqrt(power((ppt{ii}.x.position(timeOff)+ ppt{ii}.xoffset * ppt{ii}.pixelAngle),2)...
        + power((ppt{ii}.y.position(timeOff)+ ppt{ii}.yoffset * ppt{ii}.pixelAngle),2));
    if (landDist>60)
%     if (landDist>30)
        valid.landDist(ii) = true;
        counter.landDist = counter.landDist+1;   
        valid.name(ii) = 'landDist';
        valid.discard(ii) = true;
        counter.Discard = counter.Discard+1;
        continue;
    end
        
    %Drift trials
    [~,emat] = (min(abs(ppt{ii}.saccades.start - SaccON)));
    if ~isempty(emat)
        EMAT_saccOn = ppt{ii}.saccades.start(emat);
        EMAT_saccOff = EMAT_saccOn + ppt{ii}.saccades.duration(emat);
    end
    [~,emat] = (min(abs(ppt{ii}.drifts.start - SaccOFF)));
    if ~isempty(emat)
        EMAT_driftOn = ppt{ii}.drifts.start(emat);
        EMAT_driftOff = EMAT_driftOn + ppt{ii}.drifts.duration(emat);
    end
    if (EMAT_driftOn==EMAT_saccOff && EMAT_driftOff>=StimOFF && EMAT_driftOn<StimOFF)
            valid.drift(ii) = true ;
            counter.Drifts = counter.Drifts+1;
            valid.name(ii) = 'Drift';
            valid.validTrials(ii) = true ;
            counter.validTrials = counter.validTrials+1;
            continue;
    end
    
    [Answer, Intersected]=isIntersectedIn(timeOn,StimOFF-timeOn,ppt{ii}.microsaccades);
    if Answer%%%MICROSACCADE TRIALS
        if length(Intersected)==1
            valid.micsac(ii) = true;
            valid.micsac1(ii) = true;
            counter.Microsaccades = counter.Microsaccades+1;%Overall counter
            counter.Microsaccades1 = counter.Microsaccades1+1;%Specific counter
            valid.name(ii) = 'MicroSaccade1';
            valid.validTrials(ii) = true ;
            counter.validTrials = counter.validTrials+1;
            continue;
        elseif length(Intersected)==2
            valid.micsac(ii) = true;
            valid.micsac2(ii) = true;
            counter.Microsaccades = counter.Microsaccades+1;%Overall counter
            counter.Microsaccades2 = counter.Microsaccades2+1;%Specific counter
            valid.name(ii) = 'MicroSaccade2';            
            valid.validTrials(ii) = true ;
            counter.validTrials = counter.validTrials+1;
            continue
        elseif length(Intersected)>=3
            valid.micsac(ii) = true;
            valid.micsacMore(ii) = true;
            counter.Microsaccades = counter.Microsaccades+1;%Overall counter
            counter.MicrosaccadesMore = counter.MicrosaccadesMore+1;%Specific counter
            valid.name(ii) = 'MicroSaccadeMore';            
            valid.validTrials(ii) = true ;
            counter.validTrials = counter.validTrials+1;
            continue
        end
    end
    
    [Answer, Intersected]=isIntersectedIn(EMAT_driftOn,StimOFF-EMAT_driftOn,ppt{ii}.saccades);
    if Answer%%%SACCADE TRIALS
        if length(Intersected)==1
            valid.sac(ii) = true;
            counter.Saccades = counter.Saccades+1;
            valid.sac1(ii) = true;
            counter.Saccades1 = counter.Saccades1+1;
            valid.name(ii) = 'Saccade1';            
            valid.validTrials(ii) = true ;
            counter.validTrials = counter.validTrials+1;
            continue;
        elseif length(Intersected)==2
            valid.sac(ii) = true;
            counter.Saccades = counter.Saccades+1;
            valid.sac2(ii) = true;
            counter.Saccades2 = counter.Saccades2+1;
            valid.name(ii) = 'Saccade2';            
            valid.validTrials(ii) = true ;
            counter.validTrials = counter.validTrials+1;
            continue;
        elseif length(Intersected)>=3
            valid.sac(ii) = true;
            counter.Saccades = counter.Saccades+1;
            valid.sacMore(ii) = true;
            counter.SaccadesMore = counter.SaccadesMore+1;
            valid.name(ii) = 'SaccadeMore';            
            valid.validTrials(ii) = true ;
            counter.validTrials = counter.validTrials+1;
            continue;
        end
    end
    
%     if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.saccades)%%%SACCADE TRIALS
%         valid.sac(ii) = true;
%         counter.Saccades = counter.Saccades+1;
%         timeOn = round(SaccOFF);
%         timeOff = round(StimOFF);
%         valid.name(ii) = 'Saccade';
%              if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS
% %              if (isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.drifts) && ...
% %                    isIncludedInPeriod(ppt{ii}.saccades,timeOn,timeOff-timeOn))  
%         %     if isIntersectedIn(timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS
%                 valid.validTrials(ii) = true ;
%                 counter.validTrials = counter.validTrials+1;
%                 counter.Saccades = counter.Saccades-1;
%                 valid.name(ii) = 'VALID';
%                 continue;
%             end
%         continue;
%     end  
    
    %     if isIntersectedIn(timeOn,timeOff-timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS

%     [~,offline] = (min(abs(ppt{ii}.drifts.start - timeOff)));
%     if ~isempty(offline)
%         if offline<50
%             valid.drift(ii) = true ;
%             counter.Drifts = counter.Drifts+1;
%             valid.name(ii) = 'Drift';
%             continue;
%         end
%     end
%     if isIncludedInPeriod(ppt{ii}.drifts,timeOn,timeOff-timeOn)%%%DRIFT TRIALS
%         %     if isIntersectedIn(timeOn,ppt{ii}.drifts)%%%DRIFT TRIALS
%         valid.drift(ii) = true ;
%         counter.Drifts = counter.Drifts+1;
%         valid.name(ii) = 'Drift';
%         continue;
%     end

    valid.noTags(ii) = true;
    counter.noTags = counter.noTags+1;
    valid.name(ii) = 'noTags';
    valid.discard(ii) = true;
    counter.Discard = counter.Discard+1;
   
end

fprintf('\nThere are %i in total trials.\n', numTotalTrials)
fprintf('There are %i trials that have blinks.\n',counter.Blinks)
fprintf('There are %i trials do not track.\n',counter.NoTracks)
fprintf('There are %i trials that have been mis-tagged.\n',counter.misTags)
fprintf('There are %i trials that have exceeded landing distance of 60 arcmin.\n',counter.landDist)
fprintf('There are %i trials that have no tagged events\n',counter.noTags)
fprintf('There are %i trials that have been discarded.\n',counter.Discard)
fprintf('There are %i trials that have saccades in the stimulus presentation period of which:\n',counter.Saccades)
fprintf('%i trials have 1 saccade in addition to the initial saccade\n',counter.Saccades1)
fprintf('%i trials have 2 saccades in addition to the initial saccade\n',counter.Saccades2)
fprintf('%i trials have 3 or more saccades in addition to the initial saccade\n',counter.SaccadesMore)
fprintf('There are %i trials that have microsaccades in the stimulus presentation period of which:\n',counter.Microsaccades)
fprintf('%i trials have 1 MS\n',counter.Microsaccades1)
fprintf('%i trials have 2 MS\n',counter.Microsaccades2)
fprintf('%i trials have 3 or more MS\n',counter.MicrosaccadesMore)
fprintf('There are %i trials that have only drift in the post-saccade period.\n',counter.Drifts)
fprintf('There are %i trials that are valid.\n',counter.validTrials)
% fprintf('There are %i trials after filter.\n',  sum(cell2mat(struct2cell(counter))) - counter.TotalTask)