pathtodata = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data\Nikunj';
bin = 0; %set to 1 if binning
read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

[noFiles,~] = size(ppt);
trials = struct;
for i=1:noFiles
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
    
    if trials(i).contrast > 0.5
        trials(i).contrast = 0.5;
    end
    trials(i).blinks = ppt{i,1}.blinks;
    trials(i).notracks = ppt{i,1}.notracks;
    trials(i).invalid = ppt{i,1}.invalid;
    trials(i).drifts = ppt{i,1}.drifts;
    trials(i).saccades = ppt{i,1}.saccades;
    trials(i).microsaccades = ppt{i,1}.microsaccades;
end

ecc_levels = unique([trials.ecc]);
pres_levels= unique([trials.presTime]);
spFreq_lev = unique([trials.spFreq]);
% contrast = [trials.contrast];
% pres_Time = [trials.presTime];
% sp_Freq = [trials.spFreq];
conditions = struct([]);
count = 1;
for i=1:length(ecc_levels)
    for j=1:length(pres_levels)
        for k=1:length(spFreq_lev)
            In_logi =([trials.presTime]==pres_levels(j) & [trials.ecc]==ecc_levels(i)...
                & [trials.spFreq]==spFreq_lev(k) & [trials.present]==1); %indices
            In = find(In_logi);
            [valid,counter] = countingTrialsNK_mat(trials(In));
            conditions(count).ecc = ecc_levels(i);
            conditions(count).presTime = pres_levels(j);
            conditions(count).spFreq= spFreq_lev(k);
            conditions(count).validTrials = counter.validTrials;
            conditions(count).Blinks = counter.Blinks;
            conditions(count).NoTracks = counter.NoTracks;
            conditions(count).Saccades = counter.Saccades;
            conditions(count).Microsaccades = counter.Microsaccades;
            conditions(count).Drifts = counter.Drifts;
            conditions(count).manualDiscard = counter.manualDiscard;
            conditions(count).noTags = counter.noTags;
            conditions(count).totalTrials = counter.totalTrials;
            conditions(count).MS_1 = 0;
            conditions(count).MS_2 = 0;
            conditions(count).MS_more = 0;
            conditions(count).exposure = 0;
            for val=1:counter.totalTrials
                if valid.validTrials(val)
                    timeOn = round(trials(In(val)).saccOn);%Stimulus exposure period
                    timeOff = round(trials(In(val)).stimOff);
                    conditions(count).exposure = [conditions(count).exposure; timeOff-timeOn]; 
                    
                    [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).microsaccades);
                    
                    if length(Intersected)==1
                        conditions(count).MS_1 = conditions(count).MS_1 + 1;
                    elseif length(Intersected)==2
                        conditions(count).MS_2 = conditions(count).MS_2 + 1;
                    elseif length(Intersected)>2
                        conditions(count).MS_more = conditions(count).MS_more + 1;
                    end
                    
                    
%                     if Answer
%                         conditions(count).MSamp =  trials(In(val)).microsaccades.amplitude(Intersected);
%                         conditions(count).MSangle =  trials(In(val)).microsaccades.angle(Intersected);
%                     end
                    
                    timeOn = round(trials(In(val)).saccOn);%Saccade period
                    timeOff = round(trials(In(val)).saccOff);
                    [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).saccades);
                    if Answer
                        S(count_S).amp =  trials(In(val)).saccades.amplitude(Intersected);
                        S(count_S).angle =  trials(In(val)).saccades.angle(Intersected);
                        S(count_S).dur =  timeOff-timeOn;
                        S(count_S).onSpeed = trials(In(val)).velocity(timeOn);
                        %             S(count_S).vel = trials(In(val)).velocity(timeOn:timeOff);
                        S(count_S).dur_3 = sum(trials(In(val)).velocity(timeOn:timeOff) < 180);
                        S(count_S).landDist = sqrt(power(trials(In(val)).x.position(timeOff),2) + ...
                            power(trials(In(val)).y.position(timeOff),2));
                        count_S = count_S+1;
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
            count = count+1;
        end
    end
end

[valid,counter] = countingTrialsNK_mat(trials);

rownames = {'Total';'Blink';'No Track';'Saccade only';'MicroSaccade only';'Drift only';'Discarded';'Valid'};
col_1 = [counter.totalTrials;counter.Blinks;counter.NoTracks;counter.Saccades;counter.Microsaccades;counter.Drifts;...
    (counter.manualDiscard+counter.weird);counter.validTrials];
val = table(col_1,'RowNames',rownames,'VariableNames',{'Trials'});
disp(val);
writetable(val,'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\trials.csv','WriteRowNames',true,'WriteVariableNames',true);

