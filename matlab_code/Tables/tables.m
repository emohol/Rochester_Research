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
    
    trials(i).velocity = ppt{i,1}.velocity;
end

ecc_levels = unique([trials.ecc]);
pres_levels= unique([trials.presTime]);
spFreq_lev = unique([trials.spFreq]);
% contrast = [trials.contrast];
% pres_Time = [trials.presTime];
% sp_Freq = [trials.spFreq];
%%
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
            conditions(count).S_1 = 0;
            conditions(count).S_2 = 0;
            conditions(count).S_more = 0;
            conditions(count).exposure = [];
            conditions(count).vel_less_3 = [];
            conditions(count).driftOnlyTrials = 0;
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
                    timeOff = round(trials(In(val)).stimOff);
                    [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).saccades);
                    if length(Intersected)==1
                        conditions(count).S_1 = conditions(count).S_1 + 1;
                    elseif length(Intersected)==2
                        conditions(count).S_2 = conditions(count).S_2 + 1;
                    elseif length(Intersected)>2
                        conditions(count).S_more = conditions(count).S_more + 1;
                    end
                    
                    conditions(count).vel_less_3 = [conditions(count).vel_less_3; ...
                        sum(trials(In(val)).velocity(timeOn:timeOff) < 180)];%Velocity less than 3deg/s during stimulus exposure
                    
                    %                     if Answer
                    %                         S(count_S).amp =  trials(In(val)).saccades.amplitude(Intersected);
                    %                         S(count_S).angle =  trials(In(val)).saccades.angle(Intersected);
                    %                         S(count_S).dur =  timeOff-timeOn;
                    %                         S(count_S).onSpeed = trials(In(val)).velocity(timeOn);
                    %                         %             S(count_S).vel = trials(In(val)).velocity(timeOn:timeOff);
                    %                         S(count_S).dur_3 = sum(trials(In(val)).velocity(timeOn:timeOff) < 180);
                    %                         S(count_S).landDist = sqrt(power(trials(In(val)).x.position(timeOff),2) + ...
                    %                             power(trials(In(val)).y.position(timeOff),2));
                    %                         count_S = count_S+1;
                    %                     end
                    
                    %Post-saccade analysis
                    timeOn = round(trials(In(val)).saccOff);%Post-Saccade period
                    timeOff = round(trials(In(val)).stimOff);
                    [Answer_d, Intersected_d] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).drifts);        
                    [Answer_ms, Intersected_ms] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).microsaccades); 
                    [Answer_s, Intersected_s] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).saccades); 
                    if Answer_d && ~Answer_ms && ~Answer_s
                        conditions(count).driftOnlyTrials = conditions(count).driftOnlyTrials + 1;
                    end
%                     if length(Intersected)==1
%                         conditions(count).driftOnlyTrials = conditions(count).driftOnlyTrials + 1;
%                     end
                end
            end
            conditions(count).exposure_mean = mean(conditions(count).exposure);
            conditions(count).exposure_std = std(conditions(count).exposure);
            conditions(count).vel_less_3_mean = mean(conditions(count).vel_less_3);
            conditions(count).vel_less_3_std = std(conditions(count).vel_less_3);
            count = count+1;
        end
    end
end

%Catch trials condition
for j=1:length(pres_levels)
    In_logi =([trials.presTime]==pres_levels(j) & [trials.present]==0); %indices
    In = find(In_logi);
    [valid,counter] = countingTrialsNK_mat(trials(In));
    conditions(count).ecc = [];
    conditions(count).presTime = pres_levels(j);
    conditions(count).spFreq= [];
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
    conditions(count).S_1 = 0;
    conditions(count).S_2 = 0;
    conditions(count).S_more = 0;
    conditions(count).exposure = [];
    conditions(count).vel_less_3 = [];
    conditions(count).driftOnlyTrials = 0;
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
            timeOff = round(trials(In(val)).stimOff);
            [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).saccades);
            if length(Intersected)==1
                conditions(count).S_1 = conditions(count).S_1 + 1;
            elseif length(Intersected)==2
                conditions(count).S_2 = conditions(count).S_2 + 1;
            elseif length(Intersected)>2
                conditions(count).S_more = conditions(count).S_more + 1;
            end
            
            conditions(count).vel_less_3 = [conditions(count).vel_less_3; ...
                sum(trials(In(val)).velocity(timeOn:timeOff) < 180)];%Velocity less than 3deg/s during stimulus exposure
            
            %                     if Answer
            %                         S(count_S).amp =  trials(In(val)).saccades.amplitude(Intersected);
            %                         S(count_S).angle =  trials(In(val)).saccades.angle(Intersected);
            %                         S(count_S).dur =  timeOff-timeOn;
            %                         S(count_S).onSpeed = trials(In(val)).velocity(timeOn);
            %                         %             S(count_S).vel = trials(In(val)).velocity(timeOn:timeOff);
            %                         S(count_S).dur_3 = sum(trials(In(val)).velocity(timeOn:timeOff) < 180);
            %                         S(count_S).landDist = sqrt(power(trials(In(val)).x.position(timeOff),2) + ...
            %                             power(trials(In(val)).y.position(timeOff),2));
            %                         count_S = count_S+1;
            %                     end
            
            %Post-saccade analysis
            timeOn = round(trials(In(val)).saccOff);%Post-Saccade period
            timeOff = round(trials(In(val)).stimOff);
            [Answer_d, Intersected_d] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).drifts);
            [Answer_ms, Intersected_ms] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).microsaccades);
            [Answer_s, Intersected_s] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).saccades);
            if Answer_d && ~Answer_ms && ~Answer_s
                conditions(count).driftOnlyTrials = conditions(count).driftOnlyTrials + 1;
            end
%             [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,trials(In(val)).drifts);
%             if length(Intersected)==1
%                 conditions(count).driftOnlyTrials = conditions(count).driftOnlyTrials + 1;
%             end           
            
            
        end
    end
    conditions(count).exposure_mean = mean(conditions(count).exposure);
    conditions(count).exposure_std = std(conditions(count).exposure);
    conditions(count).vel_less_3_mean = mean(conditions(count).vel_less_3);
    conditions(count).vel_less_3_std = std(conditions(count).vel_less_3);
    count = count+1;
end

%% Table
formatSpec = 'Ecc-%d, Pres-%d, SpFreq-%d';
% row = [];
con_tab = table;
con_cell = {};
for i =1:14
    name = sprintf(formatSpec,conditions(i).ecc,conditions(i).presTime,conditions(i).spFreq);
    row_1 = {name,num2str(conditions(i).totalTrials),num2str(conditions(i).validTrials), num2str(conditions(i).Drifts) ,num2str(conditions(i).Saccades), num2str(conditions(i).MS_1)...
        ,num2str(conditions(i).MS_2), num2str(conditions(i).S_2), num2str(conditions(i).S_more), [num2str(mean(conditions(i).exposure)) ' ' char(177) ' ' num2str(std(conditions(i).exposure))]...
        ,[num2str(mean(conditions(i).vel_less_3)) ' ' char(177) ' ' num2str(std(conditions(i).vel_less_3))]};
    con_cell{i,1} = row_1;
    colnames = {'Condition';'TotalTrials';'ValidTrials';'Discarded_Transac_Drift';'Discarded_No_Drift';'MS_1';...
        'MS_2';'S_2';'S_3';'Exposure';'Vel_less_3'};
    con_tab_test = cell2table(con_cell{i,1},'VariableNames',colnames);
    con_tab = [con_tab; con_tab_test];
end
writetable(con_tab,'C:\Users\Ruccilab\Box\Vis_Dynamics\Figures\conditions.csv','WriteVariableNames',true);
% con_tab = table(con_tab,'VariableNames',{colnames});
% colnames = {'Condition';'ValidTrials';'Discarded_Transac_Drift';'Discarded_No_Drift';'MS_1';...
%         'MS_2';'S_2';'S_3';'Exposure';'Vel_less_3'};
% con_tab = cell2table(con_cell);
% % con_tab = cell2table(con_cell,'VariableNames',colnames);
% disp(con_tab);
% con_tab.Properties.VariableNames = {'Condition,ValidTrials,Discarded_Transac_Drift,Discarded_No_Drift,M_1,MS_2,S_2,S_3,Exposure,Vel_less_3'};

%% To display table in figure
ar = table2array(con_tab);
f = figure;
uit = uitable(f);
% d={'Total Trials',counter.totalTrials; 'Blinks', counter.Blinks;...
%     'No Tracks',counter.NoTracks; 'Saccades only', counter.Saccades;...
%     'MicroSaccades only', counter.Microsaccades; 'Drifts only',counter.Drifts;...
%     'Discarded',(counter.manualDiscard+counter.weird); 'Valid Trials',counter.validTrials};
uit.Data = ar;
uit.ColumnName = {};   
uit.RowName = {};
uit.Position = [0 0 500 500];
uit.FontSize = 10; 

