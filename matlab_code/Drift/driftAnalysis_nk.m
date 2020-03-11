path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
subject = 'A013';
% subject = 'Nikunj';
% subject = 'A092';
% subject = 'A036';
pathtodata = fullfile(path,subject);


read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call


%%
bin = 1; %set to 1 if binning
fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures';
sub_fig_path = fullfile(fig_path,subject,'drift');
if (~isfolder(sub_fig_path))
    mkdir(sub_fig_path)
end



[valid,counter] = countingTrialsNK(ppt);

trials = struct;
count = 1;

data = struct('subject', [],...
    'trialIndex', [],...
    'duration', [],'ecc',[],'presTime',[],'spFreq',[], 'curvature', [], 'speed', [], 'varHorz', [], 'varVert', [],...
    'span', [], 'diffCoeff', [], 'Dsq', [], 'DsqSingleSeg', [],...
    'FracBM', [], 'diffCoeffBoot', [], 'FracBMBoot', [], 'dsqBoot', []);
data.subject = subject;
data(1).presTime = 50;
data(2).presTime = 50;
data(3).presTime = 500;
data(4).presTime = 500;
data(5).presTime = 50;
data(6).presTime = 50;
data(7).presTime = 500;
data(8).presTime = 500;
data(9).presTime = 50;
data(10).presTime = 50;
data(11).presTime = 500;
data(12).presTime = 500;
data(1).spFreq = 2;
data(2).spFreq = 10;
data(3).spFreq = 2;
data(4).spFreq = 10;
data(5).spFreq = 2;
data(6).spFreq = 10;
data(7).spFreq = 2;
data(8).spFreq = 10;
data(9).spFreq = 2;
data(10).spFreq = 10;
data(11).spFreq = 2;
data(12).spFreq = 10;
data(1).ecc = 0;
data(2).ecc = 0;
data(3).ecc = 0;
data(4).ecc = 0;
data(5).ecc = 4;
data(6).ecc = 4;
data(7).ecc = 4;
data(8).ecc = 4;
data(9).ecc = 8;
data(10).ecc = 8;
data(11).ecc = 8;
data(12).ecc = 8;
% drift counters
driftCount_1 = 0;
driftCount_2 = 0;
driftCount_3 = 0;
driftCount_4 = 0;
driftCount_5 = 0;
driftCount_6 = 0;
driftCount_7 = 0;
driftCount_8 = 0;
driftCount_9 = 0;
driftCount_10 = 0;
driftCount_11 = 0;
driftCount_12 = 0;
fix = struct('x', [], 'y', []);

% minDuration = 100;
smoothing = 31;
% blinkBuffer = 200;
cutseg = 0;
cutseg_front = 0;
nBoots = 10;
maxSpeed = 300;
minDuration = 10;
for i=1:length(ppt)
    if valid.drift(i)
        SaccOFF = ppt{i}.saccOff; 
        StimOFF = ppt{i}.stimOff;
        [~,emat] = (min(abs(ppt{i}.drifts.start - SaccOFF)));

        if ~isempty(emat)
            driftstart = ppt{i}.drifts.start(emat);
            driftstop = driftstart + ppt{i}.drifts.duration(emat)-1;
%             driftstop = round(StimOFF);
            if round(StimOFF)>driftstop
                cutseg = 0;
            else
                cutseg = driftstop - round(StimOFF); %cut drift
            end
        end
        if (driftstop - driftstart + 1) < minDuration
            continue; %% too short
        end
        %         driftstart = ppt{i}.drifts.start(nn);
        %         driftstop = ppt{i}.drifts.duration(nn) + driftstart - 1;
        
        % x and y positions of drift segment
        x = ppt{i}.x.position(driftstart:driftstop);
        y = ppt{i}.y.position(driftstart:driftstop);
        
        
        
        if (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 2 && ppt{i}.eccentricity == 0)
    %         continue
    %         
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            
            if (mn_speed == 0)
                continue
            end

            %save results for this drift segment
            driftCount_1 = driftCount_1 + 1;
            data(1).duration(driftCount_1) = length(x);
            data(1).trialIndex(driftCount_1) = i;

            data(1).curvature(driftCount_1) = mn_cur;
            data(1).speed(driftCount_1) = mn_speed;
            data(1).varHorz(driftCount_1) = varx;
            data(1).varVert(driftCount_1) = vary;
            data(1).span(driftCount_1) = span;
% 
%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
        elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 10&& ppt{i}.eccentricity == 0)

            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            if (mn_speed == 0)
                continue
            end
            driftCount_2 = driftCount_2 + 1;
            data(2).duration(driftCount_2) = length(x);
            data(2).trialIndex(driftCount_2) = i;
            data(2).curvature(driftCount_2) = mn_cur;
            data(2).speed(driftCount_2) = mn_speed;
            data(2).varHorz(driftCount_2) = varx;
            data(2).varVert(driftCount_2) = vary;
            data(2).span(driftCount_2) = span;

%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
        elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 2&& ppt{i}.eccentricity == 0)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_3 = driftCount_3 + 1;
            data(3).duration(driftCount_3) = length(x);
            data(3).trialIndex(driftCount_3) = i;
            data(3).curvature(driftCount_3) = mn_cur;
            data(3).speed(driftCount_3) = mn_speed;
            data(3).varHorz(driftCount_3) = varx;
            data(3).varVert(driftCount_3) = vary;
            data(3).span(driftCount_3) = span;

%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
            
        elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 10 && ppt{i}.eccentricity == 0)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_4 = driftCount_4 + 1;
            data(4).duration(driftCount_4) = length(x);
            data(4).trialIndex(driftCount_4) = i;
            data(4).curvature(driftCount_4) = mn_cur;
            data(4).speed(driftCount_4) = mn_speed;
            data(4).varHorz(driftCount_4) = varx;
            data(4).varVert(driftCount_4) = vary;
            data(4).span(driftCount_4) = span;
            
%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);


        elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 2 && ppt{i}.eccentricity == 4)
            
    %         continue
    %         
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);

            if (mn_speed == 0)
                continue
            end
            
            %save results for this drift segment
            driftCount_5 = driftCount_5 + 1;
            data(5).duration(driftCount_5) = length(x);
            data(5).trialIndex(driftCount_5) = i;

            data(5).curvature(driftCount_5) = mn_cur;
            data(5).speed(driftCount_5) = mn_speed;
            data(5).varHorz(driftCount_5) = varx;
            data(5).varVert(driftCount_5) = vary;
            data(5).span(driftCount_5) = span;
% 
%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
        elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 10&& ppt{i}.eccentricity == 4)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_6 = driftCount_6 + 1;
            data(6).duration(driftCount_6) = length(x);
            data(6).trialIndex(driftCount_6) = i;
            data(6).curvature(driftCount_6) = mn_cur;
            data(6).speed(driftCount_6) = mn_speed;
            data(6).varHorz(driftCount_6) = varx;
            data(6).varVert(driftCount_6) = vary;
            data(6).span(driftCount_6) = span;

%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
        elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 2&& ppt{i}.eccentricity == 4)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_7 = driftCount_7 + 1;
            data(7).duration(driftCount_7) = length(x);
            data(7).trialIndex(driftCount_7) = i;
            data(7).curvature(driftCount_7) = mn_cur;
            data(7).speed(driftCount_7) = mn_speed;
            data(7).varHorz(driftCount_7) = varx;
            data(7).varVert(driftCount_7) = vary;
            data(7).span(driftCount_7) = span;

%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
            
        elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 10 && ppt{i}.eccentricity == 4)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg, cutseg_front,maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_8 = driftCount_8 + 1;
            data(8).duration(driftCount_8) = length(x);
            data(8).trialIndex(driftCount_8) = i;
            data(8).curvature(driftCount_8) = mn_cur;
            data(8).speed(driftCount_8) = mn_speed;
            data(8).varHorz(driftCount_8) = varx;
            data(8).varVert(driftCount_8) = vary;
            data(8).span(driftCount_8) = span;
            
%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);


        elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 2 && ppt{i}.eccentricity == 8)
            
            %         continue
            %
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg, cutseg_front,maxSpeed);
            
            if (mn_speed == 0)
                continue
            end


            %save results for this drift segment
            driftCount_9 = driftCount_9 + 1;
            data(9).duration(driftCount_9) = length(x);
            data(9).trialIndex(driftCount_9) = i;

            data(9).curvature(driftCount_9) = mn_cur;
            data(9).speed(driftCount_9) = mn_speed;
            data(9).varHorz(driftCount_9) = varx;
            data(9).varVert(driftCount_9) = vary;
            data(9).span(driftCount_9) = span;
% 
%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
        elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 10&& ppt{i}.eccentricity == 8)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg, cutseg_front,maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_10 = driftCount_10 + 1;
            data(10).duration(driftCount_10) = length(x);
            data(10).trialIndex(driftCount_10) = i;
            data(10).curvature(driftCount_10) = mn_cur;
            data(10).speed(driftCount_10) = mn_speed;
            data(10).varHorz(driftCount_10) = varx;
            data(10).varVert(driftCount_10) = vary;
            data(10).span(driftCount_10) = span;

%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
        elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 2&& ppt{i}.eccentricity == 8)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_11 = driftCount_11 + 1;
            data(11).duration(driftCount_11) = length(x);
            data(11).trialIndex(driftCount_11) = i;
            data(11).curvature(driftCount_11) = mn_cur;
            data(11).speed(driftCount_11) = mn_speed;
            data(11).varHorz(driftCount_11) = varx;
            data(11).varVert(driftCount_11) = vary;
            data(11).span(driftCount_11) = span;

%             fix(driftCount).x = x(1+cutseg:end-cutseg);
%             fix(driftCount).y = y(1+cutseg:end-cutseg);
            
            
        elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 10 && ppt{i}.eccentricity == 8)
            
            if length(x) < smoothing
                disp('drift segment smaller than 31ms');
                continue
            end
            [span, mn_speed, mn_cur, varx, vary] = ...
                getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
            
            if (mn_speed == 0)
                continue
            end
            driftCount_12 = driftCount_12 + 1;
            data(12).duration(driftCount_12) = length(x);
            data(12).trialIndex(driftCount_12) = i;
            data(12).curvature(driftCount_12) = mn_cur;
            data(12).speed(driftCount_12) = mn_speed;
            data(12).varHorz(driftCount_12) = varx;
            data(12).varVert(driftCount_12) = vary;
            data(12).span(driftCount_12) = span;
        end
        
    elseif valid.validTrials(i)
        timeOn = ppt{i}.saccOn; 
        timeOff = ppt{i}.stimOff;
        [Answer, Intersected] = isIntersectedIn(timeOn,timeOff-timeOn,ppt{i}.drifts);
        if Answer
            for j=1:length(Intersected)
                driftstart = ppt{i}.drifts.start(j);
                driftstop = driftstart + ppt{i}.drifts.duration(j)-1;
                
                if round(timeOff)>driftstop
                    cutseg = 0;
                else
                    cutseg = driftstop - round(timeOff); %cut drift
                end
                % x and y positions of drift segment
                x = ppt{i}.x.position(driftstart:driftstop);
                y = ppt{i}.y.position(driftstart:driftstop);
                
                
                
                if (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 2 && ppt{i}.eccentricity == 0)
                    
                    %         continue
                    %
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    %save results for this drift segment
                    driftCount_1 = driftCount_1 + 1;
                    data(1).duration(driftCount_1) = length(x);
                    data(1).trialIndex(driftCount_1) = i;
                    
                    data(1).curvature(driftCount_1) = mn_cur;
                    data(1).speed(driftCount_1) = mn_speed;
                    data(1).varHorz(driftCount_1) = varx;
                    data(1).varVert(driftCount_1) = vary;
                    data(1).span(driftCount_1) = span;
                    %
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 10&& ppt{i}.eccentricity == 0)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_2 = driftCount_2 + 1;
                    data(2).duration(driftCount_2) = length(x);
                    data(2).trialIndex(driftCount_2) = i;
                    data(2).curvature(driftCount_2) = mn_cur;
                    data(2).speed(driftCount_2) = mn_speed;
                    data(2).varHorz(driftCount_2) = varx;
                    data(2).varVert(driftCount_2) = vary;
                    data(2).span(driftCount_2) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 2&& ppt{i}.eccentricity == 0)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_3 = driftCount_3 + 1;
                    data(3).duration(driftCount_3) = length(x);
                    data(3).trialIndex(driftCount_3) = i;
                    data(3).curvature(driftCount_3) = mn_cur;
                    data(3).speed(driftCount_3) = mn_speed;
                    data(3).varHorz(driftCount_3) = varx;
                    data(3).varVert(driftCount_3) = vary;
                    data(3).span(driftCount_3) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                    
                elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 10 && ppt{i}.eccentricity == 0)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_4 = driftCount_4 + 1;
                    data(4).duration(driftCount_4) = length(x);
                    data(4).trialIndex(driftCount_4) = i;
                    data(4).curvature(driftCount_4) = mn_cur;
                    data(4).speed(driftCount_4) = mn_speed;
                    data(4).varHorz(driftCount_4) = varx;
                    data(4).varVert(driftCount_4) = vary;
                    data(4).span(driftCount_4) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                    
                elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 2 && ppt{i}.eccentricity == 4)
                    
                    %         continue
                    %
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    
                    %save results for this drift segment
                    driftCount_5 = driftCount_5 + 1;
                    data(5).duration(driftCount_5) = length(x);
                    data(5).trialIndex(driftCount_5) = i;
                    
                    data(5).curvature(driftCount_5) = mn_cur;
                    data(5).speed(driftCount_5) = mn_speed;
                    data(5).varHorz(driftCount_5) = varx;
                    data(5).varVert(driftCount_5) = vary;
                    data(5).span(driftCount_5) = span;
                    %
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 10&& ppt{i}.eccentricity == 4)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_6 = driftCount_6 + 1;
                    data(6).duration(driftCount_6) = length(x);
                    data(6).trialIndex(driftCount_6) = i;
                    data(6).curvature(driftCount_6) = mn_cur;
                    data(6).speed(driftCount_6) = mn_speed;
                    data(6).varHorz(driftCount_6) = varx;
                    data(6).varVert(driftCount_6) = vary;
                    data(6).span(driftCount_6) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 2&& ppt{i}.eccentricity == 4)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_7 = driftCount_7 + 1;
                    data(7).duration(driftCount_7) = length(x);
                    data(7).trialIndex(driftCount_7) = i;
                    data(7).curvature(driftCount_7) = mn_cur;
                    data(7).speed(driftCount_7) = mn_speed;
                    data(7).varHorz(driftCount_7) = varx;
                    data(7).varVert(driftCount_7) = vary;
                    data(7).span(driftCount_7) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                    
                elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 10 && ppt{i}.eccentricity == 4)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_8 = driftCount_8 + 1;
                    data(8).duration(driftCount_8) = length(x);
                    data(8).trialIndex(driftCount_8) = i;
                    data(8).curvature(driftCount_8) = mn_cur;
                    data(8).speed(driftCount_8) = mn_speed;
                    data(8).varHorz(driftCount_8) = varx;
                    data(8).varVert(driftCount_8) = vary;
                    data(8).span(driftCount_8) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                    
                elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 2 && ppt{i}.eccentricity == 8)
                    
                    %         continue
                    %
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    
                    %save results for this drift segment
                    driftCount_9 = driftCount_9 + 1;
                    data(9).duration(driftCount_9) = length(x);
                    data(9).trialIndex(driftCount_9) = i;
                    
                    data(9).curvature(driftCount_9) = mn_cur;
                    data(9).speed(driftCount_9) = mn_speed;
                    data(9).varHorz(driftCount_9) = varx;
                    data(9).varVert(driftCount_9) = vary;
                    data(9).span(driftCount_9) = span;
                    %
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                elseif (ppt{i}.presTime == 50 && ppt{i}.spatialFreq == 10&& ppt{i}.eccentricity == 8)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_10 = driftCount_10 + 1;
                    data(10).duration(driftCount_10) = length(x);
                    data(10).trialIndex(driftCount_10) = i;
                    data(10).curvature(driftCount_10) = mn_cur;
                    data(10).speed(driftCount_10) = mn_speed;
                    data(10).varHorz(driftCount_10) = varx;
                    data(10).varVert(driftCount_10) = vary;
                    data(10).span(driftCount_10) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 2&& ppt{i}.eccentricity == 8)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_11 = driftCount_11 + 1;
                    data(11).duration(driftCount_11) = length(x);
                    data(11).trialIndex(driftCount_11) = i;
                    data(11).curvature(driftCount_11) = mn_cur;
                    data(11).speed(driftCount_11) = mn_speed;
                    data(11).varHorz(driftCount_11) = varx;
                    data(11).varVert(driftCount_11) = vary;
                    data(11).span(driftCount_11) = span;
                    
                    %             fix(driftCount).x = x(1+cutseg:end-cutseg);
                    %             fix(driftCount).y = y(1+cutseg:end-cutseg);
                    
                    
                elseif (ppt{i}.presTime == 500 && ppt{i}.spatialFreq == 10 && ppt{i}.eccentricity == 8)
                    
                    if length(x) < smoothing
                        disp('drift segment smaller than 31ms');
                        continue
                    end
                    [span, mn_speed, mn_cur, varx, vary] = ...
                        getDriftChar_nk(x, y, smoothing, cutseg,cutseg_front, maxSpeed);
                    
                    if (mn_speed == 0)
                        continue
                    end
                    driftCount_12 = driftCount_12 + 1;
                    data(12).duration(driftCount_12) = length(x);
                    data(12).trialIndex(driftCount_12) = i;
                    data(12).curvature(driftCount_12) = mn_cur;
                    data(12).speed(driftCount_12) = mn_speed;
                    data(12).varHorz(driftCount_12) = varx;
                    data(12).varVert(driftCount_12) = vary;
                    data(12).span(driftCount_12) = span;
                end
                
            end
        end
    end
end


%% Plot
fh = 0;
for i=1:12
    allDSpeed = [];
    allDSpeed = cat(2, data(i).speed);
    ntotal = length(allDSpeed);
    nbins = max(10, min(60, round(ntotal/3)));
    % c = linspace(0, 60, nbins);
    c = linspace(min(allDSpeed), max(allDSpeed), nbins);
    m1 = nanmean(allDSpeed(:));
    s1 = nanstd(allDSpeed(:));% / sqrt(sum(~isnan(allMSAmps(:))));
    n1 = hist(allDSpeed(:), c);
    n1 = n1 / sum(n1);    
    fh=fh+1;
    figure(fh)
    % hist(allDSpeed(:), c);
    histogram(allDSpeed(:), c);
    % plot(c, n1,'linewidth', 2);
    vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
    set(gca,'box','off','tickdir','out','FontSize',20)
    xlabel('Speed (arcmin/sec)');
    ylabel('Prob');
    % xlim([0 10])
    formatSpec = "Mean-%1.1f,Std-%1.1f";
    txt = sprintf(formatSpec,m1,s1);
%     text(round(m1),round(max(n1)),txt,'FontSize',20);
    text(double(round(m1-s1)),0,txt,'FontSize',20);
    saveas(figure(fh),[sub_fig_path '\D-Speed'],'epsc');
    saveas(figure(fh),[sub_fig_path '\Ecc',num2str(data(i).ecc),'_Pres_',num2str(data(i).presTime),'_SpFreq_',num2str(data(i).spFreq),'_D-Speed'],'epsc');
    
    allDCur = [];
    allDCur = cat(2, data(i).curvature);
    ntotal = length(allDCur);
    nbins = max(10, min(60, round(ntotal/3)));
    % c = linspace(0, 60, nbins);
    c = linspace(min(allDCur), max(allDCur), nbins);
    m1 = nanmean(allDCur(:));
    s1 = nanstd(allDCur(:));% / sqrt(sum(~isnan(allMSAmps(:))));
    n1 = hist(allDCur(:), c);
    n1 = n1 / sum(n1);
    fh = fh+1;
    figure(fh)
    histogram(allDCur(:), c);
    % plot(c, n1,'linewidth', 2);
    vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
    set(gca,'box','off','tickdir','out','FontSize',20)
    xlabel('Curvature');
    ylabel('Prob');
    xlim([-50 500]);
    % xlim([min(data.curvature) max(data.curvature)])
    formatSpec = "Mean-%1.1f,Std-%1.1f";
    txt = sprintf(formatSpec,m1,s1);
    text(double(round(m1)),0,txt,'FontSize',20);
%     text(m1,max(n1),txt,'FontSize',20);
   saveas(figure(fh),[sub_fig_path '\Ecc',num2str(data(i).ecc),'_Pres_',num2str(data(i).presTime),'_SpFreq_',num2str(data(i).spFreq),'_D-Cur'],'epsc');
    
   allDSpan = [];
    allDSpan = cat(2, data(i).span);
    ntotal = length(allDSpan);
    nbins = max(10, min(60, round(ntotal/3)));
    % c = linspace(0, 60, nbins);
    c = linspace(min(allDSpan), max(allDSpan), nbins);
    m1 = nanmean(allDSpan(:));
    s1 = nanstd(allDSpan(:));% / sqrt(sum(~isnan(allMSAmps(:))));
    n1 = hist(allDSpan(:), c);
    n1 = n1 / sum(n1);
    fh = fh+1;
    figure(fh)
    histogram(allDSpan(:), c);
    vertLineThrough([m1, m1+s1, m1-s1],'r', gca ,'--');
    set(gca,'box','off','tickdir','out','FontSize',20)
    xlabel('Span');
    ylabel('Prob');
    formatSpec = "Mean-%1.1f,Std-%1.1f";
    txt = sprintf(formatSpec,m1,s1);
    text(double(round(m1-s1)),0,txt,'FontSize',20);
%     text(m1,max(n1),txt,'FontSize',20);
    saveas(figure(fh),[sub_fig_path '\Ecc',num2str(data(i).ecc),'_Pres_',num2str(data(i).presTime),'_SpFreq_',num2str(data(i).spFreq),'_D-Span'],'epsc');
end
%%
% diffusion coefficient estimation and bootstrapping
[~, ~, dc, ~, dsq, dsqss] = CalculateDiffusionCoef(fix, 'MaxTime', max(data.duration));
%fitpar = fit(log(4*(50:255)/1000)', dsq(50:end)', 'exp1');
fitpar = struct('a', 0, 'b', 0);

data.diffCoeff = dc;
data.Dsq = dsq;
data.DsqSingleSeg = dsqss;
data.FracBM = [fitpar.a, fitpar.b];

tempdiffCoeffBoot = nan(nBoots, 1);
tempFracBMBoot = nan(nBoots, 2);
tempDsq = nan(nBoots, length(dsq));
dsqnn = length(dsq);
parfor nb = 1:nBoots
    tmp = randsample(length(fix), length(fix), true);
    fix2 = fix(tmp);
    [~, ~, dc, ~, dsqcc, ~] = CalculateDiffusionCoef(fix2, 'MaxTime', max(data.duration));
    %fitpar = fit(log(4*(50:255)/1000)', dsq(50:end)', 'exp1');
    fitpar = struct('a', 0, 'b', 0);
    
    tempdiffCoeffBoot(nb) = dc;
    tempFracBMBoot(nb, :) = [fitpar.a, fitpar.b];
    
    ll = min(length(dsqcc), dsqnn);
    dsqtt = nan(1, dsqnn);
    dsqtt(1:ll) = dsqcc(1:ll);
    tempDsq(nb, :) = dsqtt;
end

data.diffCoeffBoot = tempdiffCoeffBoot;
data.FracBMBoot = tempFracBMBoot;
data.dsqBoot = tempDsq;