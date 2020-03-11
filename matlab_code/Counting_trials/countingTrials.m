% function [counter, valid] = countingTrials(vt, trials)%,numTooShort,numBlinks,numNoTracks,numSaccades,numMicroSaccades)
%COUNTING TRIALS Summary of this function goes here
%   Detailed explanation goes here
pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj_NoMask';
% pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj';
read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[pptrials] = preprocessing(data); % this is a standard preprocessing call
numTotalTrials = length(pptrials);
% dist_thresh = 30;

% pptrials = osRemoverInEM(pptrials); %To use or not?

counter.TooShort = 0;
counter.BigOffset = 0;
counter.Blinks = 0;
counter.NoTracks = 0;
counter.NoResponse = 0;
counter.manualDiscard = 0;
counter.Saccades = 0;
counter.Fixation = 0;
counter.Span = 0;
counter.Microsaccades = 0;
counter.Drifts = 0;
counter.TotalTask = 0;

valid = struct(...;
    'tooshort', false(size(pptrials)),...;
    'blink', false(size(pptrials)),...;
    'notrack', false(size(pptrials)),...;
    'noresponse', false(size(pptrials)),...;
    'manualdiscard', false(size(pptrials)),...;
    's', false(size(pptrials)),...;
    'fixation', false(size(pptrials)),...;
    'span', false(size(pptrials)),...;
    'dms', false(size(pptrials)),...;
    'ms', false(size(pptrials)),...;
    'drift', false(size(pptrials)));

%different counts for different events
for ii = 1:length(pptrials)
%     if length(pptrials{ii}.x.position) < 500 %%% TOO SHORT
%         valid.tooshort(ii) = true;
%         counter.TooShort = counter.TooShort+1;
%         continue;
%     end
%     if pptrials{ii}.Correct < 3 % TASK TRIAL
%         if figures.FIXATION_ANALYSIS %%%Pre-Target Time ~15ms
%             timeOn = round(pptrials{ii}.TimeCueON);
%             timeOff = round(pptrials{ii}.TimeCueOFF);
%             x = pptrials{ii}.x.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.xoffset;
%             y = pptrials{ii}.y.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.yoffset;
%             em.span(ii) = quantile(sqrt((x-mean(x)).^2 + (y-mean(y)).^2), .95);
%             
%         else %%%When Target is presented ~500ms
%             
%             timeOn = round(pptrials{ii}.TimeTargetON);
%             timeOff = round(min(pptrials{ii}.TimeTargetOFF, pptrials{ii}.ResponseTime));
%             x = pptrials{ii}.x.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.xoffset;
%             y = pptrials{ii}.y.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.yoffset;
%             em.span(ii) = quantile(sqrt((x-mean(x)).^2 + (y-mean(y)).^2), .95);
%         end
%         
%     else % FIXATION TRIAL
%         timeOn = round(pptrials{ii}.TimeFixationON)+1;
%         timeOff = round(pptrials{ii}.TimeTargetOFF);
%         x = pptrials{ii}.x.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.xoffset;
%         y = pptrials{ii}.y.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.yoffset;
%         em.span(ii) = quantile(sqrt((x-mean(x)).^2 + (y-mean(y)).^2), .95);
%     end
%     
%     if sum(abs(pptrials{ii}.x.position(timeOn:(timeOff-timeOn))... %%%BIG-OFFSET 
%          + pptrials{ii}.xoffset * pptrials{ii}.pxAngle) > dist_thresh) > 0 || ...
%        sum(abs(pptrials{ii}.y.position(timeOn:(timeOff-timeOn))...
%          + pptrials{ii}.yoffset * pptrials{ii}.pxAngle) > dist_thresh) > 0
%             pptrials{ii}.BigOffset = 1;
%             valid.BigOffset(ii) = true;
%             counter.BigOffset = counter.BigOffset+1;
%         continue;
%     end

%     timeOn = round(pptrials{ii}.TimeCueON);
%     timeOff = round(pptrials{ii}.TimeCueOFF);
%     x = pptrials{ii}.x.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.xoffset;
%     y = pptrials{ii}.y.position(timeOn:timeOff) + params.pixelAngle * pptrials{ii}.yoffset;
    CuePrf = pptrials{ii}.cuePrf;
    FixationON = pptrials{ii}.fixOn;
    CueCtr = pptrials{ii}.cueCtr;
%     FixationOFF = FixationON + pptrials{ii}.fixTime;     
    SaccON = pptrials{ii}.saccOn;
    FlashON = pptrials{ii}.flashOn;
    SaccOFF = pptrials{ii}.saccOff;
%     TargetON = pptrials{ii}.TimeTargetON;    
    StimOFF = pptrials{ii}.stimOff;
    Quit = pptrials{ii}.quit;
    RespTime = pptrials{ii}.responseTime;
%     RespCueOFF = RespCueON + pptrials{ii}.CueTime; 
%     ResponseTime = pptrials{ii}.ResponseTime;
%     timeOn = round(FlashON);
    timeOn = round(SaccOFF);%Cause presentation time starts after central fixation
%     timeOff = round(RespTime);
    timeOff = round(StimOFF);
    if isIntersectedIn(timeOn,timeOff-timeOn,pptrials{ii}.blinks)%%%BLINKS
        valid.blink(ii) = true;
        counter.Blinks = counter.Blinks+1;
        continue;
    end
    
    if isIntersectedIn(timeOn,timeOff-timeOn,pptrials{ii}.notracks)%%%NO TRACK TRIALS
        valid.notrack(ii) = true;
        counter.NoTracks = counter.NoTracks+1;
        continue;
    end
%     if pptrials{ii}.Correct == 3 %%%NO RESPONSE
%        valid.noresponse(ii) = true;
%        counter.NoResponse = counter.NoResponse+1;
%        continue;
%     end
    if  ~isempty(pptrials{ii}.invalid.start)
        valid.manualdiscard(ii) = true;
        counter.manualDiscard = counter.manualDiscard+1;
        continue;
    end
    if isIntersectedIn(timeOn,timeOff-timeOn,pptrials{ii}.saccades)%%%SACCADE TRIALS
        valid.s(ii) = true;
        counter.Saccades = counter.Saccades+1;
        continue;
    end
%     if pptrials{ii}.TimeFixationOFF > 0 %%%FIXATION TRIALS
%         valid.fixation(ii) = true;
%         counter.Fixation = counter.Fixation+1;
%         continue;
%     end
%     if ~(params.spanMin < em.span(ii) && em.span(ii) < params.spanMax)
%         valid.span(ii) = true;
%         counter.Span = counter.Span+1;
%         continue;
%     end    
    if isIntersectedIn(timeOn,timeOff-timeOn,pptrials{ii}.microsaccades)%%%MICROSACCADE TRIALS
        valid.ms(ii) = true;
        valid.dms(ii) = true;
        counter.Microsaccades = counter.Microsaccades+1;
        counter.TotalTask = counter.TotalTask + 1;
        continue;
%     else 
%         valid.drift(ii) = true;
%         valid.dms(ii) = true;
%         counter.Drifts = counter.Drifts+1;
%         counter.TotalTask = counter.TotalTask + 1;
    end
    if isIntersectedIn(timeOn,timeOff-timeOn,pptrials{ii}.drifts)%%%DRIFT TRIALS
%     if isIntersectedIn(timeOn,pptrials{ii}.drifts)%%%DRIFT TRIALS
        valid.drift(ii) = true ;
        counter.Drifts = counter.Drifts+1;
        continue;
    end    
end

fprintf('\nThere are %i in total trials prior to filter.\n\n', numTotalTrials)
fprintf('There are %i trials that have no response.\n',counter.TooShort)
fprintf('There are %i trials that have large offset.\n',counter.BigOffset)
fprintf('There are %i trials that have blinks.\n',counter.Blinks)
fprintf('There are %i trials do not track.\n',counter.NoTracks)
fprintf('There are %i trials that have no response.\n',counter.NoResponse)
fprintf('There are %i trials that have saccades.\n',counter.Saccades)
fprintf('There are %i trials that are fixation.\n',counter.Fixation)
fprintf('There are %i trials that do not fit within the span.\n',counter.Span)
fprintf('There are %i trials that have microsaccades.\n',counter.Microsaccades)
fprintf('There are %i trials that have drifts only.\n',counter.Drifts)
fprintf('There are %i trials after filter.\n',  sum(cell2mat(struct2cell(counter))) - counter.TotalTask)

% end



