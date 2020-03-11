%% notes to self
% monitor parameters (ASUS 278 at 160cm away)
%       1440 x 900 pixels @ 144 Hz (~0.8 arcmin per pixel)
%       asus brightness and contrast = 0
%       nvidia gamma corrections: 2.1, 2.17, 2.58 (RGB)

%%
% pptrials is a cell array of structs with the following fields:
% (Check the Data/subject folders in the google drive for pptrials.mat)
%
%   *** General Experiment Data ***
%       Subject_Name - string
%       stimulusEnvelope - should be gaussian.tga
%
%   *** Block Parameters ***
%       spatialFreq - 6 or 20 cpd
%       targetPerf - 0.65 or 0.9 proportions
%
%   *** Trial Parameters ***
%       contrast - float (grating with this amplitude)
%       orientation - 0 or 1
%
%   *** Subject Response Data ***
%       Response - from 0 or 1 (1000 if no response)
%       Correct - 0 or 1, (3 if no response)
%       ResponseTime - time (ms) in trial of button press
%
%   *** Stimulus Presentation Timing (all in ms) ***
%       TimeHoldON, TimeHoldOFF      (blank period before stimulus)
%       TimeRampON, TimeRampOFF      (stimulus contrast ramp)
%       TimePlatON, TimePlatOFF      (stimulus contrast plateau)
%
%   Other fields are related to fixation/calibration trials or contain
%   eye movement data.

%% run parameters
subj = 'A047_DDPI'; % A031, A039, A084, A047_DDPI

imgPath = fullfile(pathToBOX(), 'DriftControl', 'Figures', subj);
if ~exist(imgPath, 'dir'), mkdir(imgPath), end

% specify all parameters in ms - need to convert to samples elsewhere
diffConDuration = 512;
nBoots = 2;
nContrastGroups = 3;

pixelAngle = .8;

gray2MichContrast = @(gray) gray / 127;

switch subj
    case 'A031'
        analysisDuration = 512;
    case {'A084', 'A039'}
        analysisDuration = 1024;
    otherwise % new subject
        analysisDuration = 512;
end

%% plotting parameters
align = {'right', 'left'};
thrCols = 'gm';
spFCols = 'br';

%% load in pptrials
pt = fullfile(pathToOPUS(), 'DriftControl_spf', subj, 'pptrials.mat');

spFreq = [6, 20]; % list of spatial frequencies to analyze
targetPerf = [.61, .89]; % target performances to analyze
perfWidth = .12;
% load subject data
load(fullfile(pt), 'pptrials', 'data');
subj = pptrials{1}.Subject_Name;
counts = struct(...
    'notracks', false(length(pptrials), 1),...
    'blinks', false(length(pptrials), 1),...
    'saccades', false(length(pptrials), 1),...
    'microsaccades', false(length(pptrials), 1),...
    'driftonly', false(length(pptrials), 1),...
    'fixation', false(length(pptrials), 1));

%% initialize directories for saving things
savePt = fileparts(pt); % just get the directory of the path - save stuff here
if ~exist(fullfile(savePt, 'CountTables'), 'dir')
    mkdir( fullfile(savePt, 'CountTables') );
end

fh = 0; % figure handles
%% rearrange pptrials into different data struct
flds = fields(data.user{1});
flds(4) = []; % remove Subject_Name
flds(end-3) = []; % remove stimulusEnvelope

clear data;
data = struct();

for fi = 1:length(flds)
    eval(sprintf('data.%s = cellfun(@(x) x.%s, pptrials);', flds{fi}, flds{fi}));
end

% round contrast values to 2 decimal points - this improves the binning for
% the psychometric function fitting
data.contrast = round(data.contrast, 2); % raw contrast
data.michContrast = gray2MichContrast(data.contrast); % convert to michelson contrast
data.michContrastRounded = 10.^round(log10(data.michContrast), 1);
data.spFreq = spFreq;
data.targetPerf = targetPerf;

%% count valid trials
for ii = 1:length(pptrials) % loop through each trial
    
    if pptrials{ii}.TimeFixationON == 0 % these are in time, convert to samples
        stimStop = round(pptrials{ii}.TimePlatOFF) / 1000 * pptrials{ii}.Fs;
        stimStart = stimStop - analysisDuration / 1000 * pptrials{ii}.Fs;
    else
        stimStart = pptrials{ii}.TimeFixationON / 1000 * pptrials{ii}.Fs;
        stimStop = pptrials{ii}.TimeFixationOFF / 1000 * pptrials{ii}.Fs;
        
        counts.fixation(ii) = true;
    end
    
    % some trials just don't get classified - check these out:
    % A084 - 250, 348
    if isIntersectedIn(stimStart, stimStop - stimStart, pptrials{ii}.notracks)|| ...
            isempty(pptrials{ii}.drifts.start) || ... % bad tracking
            manualFilter(subj, ii)
        counts.notracks(ii) = true;
    elseif isIntersectedIn(stimStart, stimStop - stimStart, pptrials{ii}.blinks)
        counts.blinks(ii) = true;
    elseif isIntersectedIn(stimStart, stimStop - stimStart, pptrials{ii}.saccades)
        counts.saccades(ii) = true;
    elseif isIntersectedIn(stimStart, stimStop - stimStart, pptrials{ii}.microsaccades)
        counts.microsaccades(ii) = true;
    else
        counts.driftonly(ii) = true;
    end
end % done looping through trials

%% loop through trials and define driftonly segments
clear fix;
fix(length(pptrials), 1) = struct('x', [], 'y', [], 'trialid', []);
cnt = 0;
for ii = 1:length(pptrials)
    x = pptrials{ii}.x.position;
    isdrift = false(size(x));
    ds = pptrials{ii}.drifts.start;
    for dd = 1:length(ds)
        sss = pptrials{ii}.drifts.start(dd);
        eee = sss + pptrials{ii}.drifts.duration(dd) - 1;
        isdrift(sss:eee) = true;
    end
    pptrials{ii}.isdrift = isdrift;
    
    %% create single fix struct  of drift only segments in diffCon period
    if pptrials{ii}.TimeFixationON == 0
        stimStop = round(pptrials{ii}.TimePlatOFF / 1000 * pptrials{ii}.Fs);
        stimStart = stimStop - round(diffConDuration / 1000 * pptrials{ii}.Fs);
    else
        stimStart = max(ceil(round(pptrials{ii}.TimeFixationON / 1000 * pptrials{ii}.Fs)), 1);
        stimStop = min(floor(round(pptrials{ii}.TimeFixationOFF / 1000 * pptrials{ii}.Fs)), length(x));
    end
    
    if counts.blinks(ii) || counts.notracks(ii)
        continue;
    end
    
    [sss, eee] = getIndicesFromBin(pptrials{ii}.isdrift(stimStart:stimStop));
    xx = pptrials{ii}.x.position(stimStart:stimStop)...
        - pptrials{ii}.xoffset*pixelAngle;
    yy = pptrials{ii}.y.position(stimStart:stimStop)...
        - pptrials{ii}.yoffset*pixelAngle;
    for kk = 1:length(sss)
        
        si = sss(kk) + 50;
        ei = eee(kk);
        
        if si >= ei
            continue;
        end
        
        cnt = cnt + 1;
        fix(cnt).x = xx(sss(kk):eee(kk));
        fix(cnt).y = yy(sss(kk):eee(kk));
        fix(cnt).trialid = ii;
    end
end
fixtrialid = [fix.trialid];
fixMichContrasts = data.michContrast(fixtrialid);
fixspf = data.spatialFreq(fixtrialid);


%% analyze drift during fixation trials
useTrials = counts.fixation & ~(counts.notracks | counts.blinks);

[~, useTrialst] = intersect(fixtrialid, find(useTrials));
fixtemp = fix(useTrialst);

[~, ~, dc, ~, dsq, dsqSingleSeg] = CalculateDiffusionCoef(fixtemp,...
    'Fsampling', pptrials{1}.Fs, 'MaxTime', 256 / 1000 * pptrials{1}.Fs);
data.fixation.diffConstant.dc = dc;
data.fixation.diffConstant.dsq = dsq;
data.fixation.diffConstant.dsqSingleSeg = dsqSingleSeg;

data.fixation.diffConstant.dcBoot = nan(nBoots, 1);
for bb = 1:nBoots
    rs = randsample(length(fixtemp), length(fixtemp), true);
    [~, ~, dc] = CalculateDiffusionCoef(fixtemp(rs),...
        'Fsampling', pptrials{1}.Fs, 'MaxTime', 256 / 1000 * pptrials{1}.Fs);
    data.fixation.diffConstant.dcBoot(bb) = dc;
end

%% microsaccade and saccade occurences / amplitudes / angles
fh = 100;
winsize = 50;
c = (-2000:winsize:1000);

microsaccadeTimes = cell(size(pptrials));
saccadeTimes = cell(size(pptrials));

microsaccadeAmps = cell(size(pptrials));
saccadeAmps = cell(size(pptrials));

microsaccadeDirs = cell(size(pptrials));
saccadeDirs = cell(size(pptrials));
for ii = 1:length(pptrials)
    if counts.blinks(ii) || counts.notracks(ii)
        continue;
    end
    if ~counts.fixation(ii)
        stimOFF = pptrials{ii}.TimePlatOFF;
    else
        stimOFF = pptrials{ii}.TimeFixationOFF;
    end
    
    % here convert to time (ms)
    microsaccadeTimes{ii} = (pptrials{ii}.microsaccades.start + ...
        pptrials{ii}.microsaccades.duration / 2) / pptrials{ii}.Fs * 1000 - ...
        stimOFF;
    saccadeTimes{ii} = (pptrials{ii}.saccades.start + ...
        pptrials{ii}.saccades.duration / 2) / pptrials{ii}.Fs * 1000 - ...
        stimOFF;
    
    microsaccadeAmps{ii} = pptrials{ii}.microsaccades.amplitude;
    microsaccadeDirs{ii} = pptrials{ii}.microsaccades.angle;
    
    saccadeAmps{ii} = pptrials{ii}.saccades.amplitude;
    saccadeDirs{ii} = pptrials{ii}.saccades.angle;
end

for ss = 1:length(spFreq)+1
    useTrials = ~(counts.blinks(:) | counts.notracks(:));
    if ss <= length(spFreq)
        useTrials = useTrials(:) & data.spatialFreq(:) == spFreq(ss);
    else
        useTrials = useTrials(:) & counts.fixation(:);
    end
    allMSTimes = cat(2, microsaccadeTimes{useTrials});
    allSTimes = cat(2, saccadeTimes{useTrials});
    allMSAmps = cat(2, microsaccadeAmps{useTrials});
    allSAmps = cat(2, saccadeAmps{useTrials});
    allMSDirs = cat(2, microsaccadeDirs{useTrials});
    allSDirs = cat(2, saccadeDirs{useTrials});
    
    allAmps = [allMSAmps, allSAmps];
    allTimes = [allMSTimes, allSTimes];
    
    %%%%%%%%% saccade/microsaccade rates %%%%%%%%%%%%%
    cnt = sum(useTrials);
    nMS = hist(allMSTimes, c);
    nS = hist(allSTimes, c);
    nMS = nMS / (winsize / 1000 * cnt);
    nS = nS / (winsize / 1000 * cnt);
    
    fh = fh + 1;
    figure(fh); clf; hold on;
    hhh(1) = plot(c(2:end-1), nMS(2:end-1), 'g', 'linewidth', 2);
    hhh(2) = plot(c(2:end-1), nS(2:end-1), 'b', 'linewidth', 2);
    hhh(3) = plot(c(2:end-1), nMS(2:end-1) + nS(2:end-1), 'k', 'linewidth', 2);
    vertLineThrough(0, 'k', gca, ':');
    
    xlabel('time (ms) relative to stimulus offset');
    ylabel('rate (# saccade / s)');
    legend(hhh, {'microsaccade', 'saccade', 'all'}, 'Location', 'east');
    ylim([0, 3]);
    yl = ylim;
    
    if ss <= length(spFreq)
        vertLineThrough(-1300, 'k', gca, ':');
        vertLineThrough(-500, 'k', gca, ':');
        text(-1000, yl(1) + .9 * diff(yl), 'Ramp',...
            'FontSize', 14, 'FontWeight', 'bold',...
            'HorizontalAlignment', 'center');
        text(-300, yl(1) + .9 * diff(yl), 'Plat',...
            'FontSize', 14, 'FontWeight', 'bold',...
            'HorizontalAlignment', 'center');
        text(500, yl(1) + .9 * diff(yl), 'Mask',...
            'FontSize', 14, 'FontWeight', 'bold',...
            'HorizontalAlignment', 'center');
    end
    
    standardPlot;
    if ss <= length(spFreq)
        title(sprintf('%s - %icpd', strrep(subj, '_', ' '), spFreq(ss)));
        print(fullfile(imgPath, sprintf('%s_saccadeRate_%icpd.png', subj, spFreq(ss))), '-dpng');
        print(fullfile(imgPath, sprintf('%s_saccadeRate_%icpd.eps', subj, spFreq(ss))), '-depsc');
    else
        title(sprintf('%s - Fixation', strrep(subj, '_', ' '), subj));
        print(fullfile(imgPath, sprintf('%s_saccadeRate_Fixation.png', subj)), '-dpng');
        print(fullfile(imgPath, sprintf('%s_saccadeRate_Fixation.eps', subj)), '-depsc');
    end
    
end

%% initialize structures for saving data
allowList = {'driftonly', 'microsaccades', 'saccades', 'ms_only', 'ms_saccades_only'};

data.targetPerf = targetPerf;

% initialize structures for psychometric fitting info
data.psyfit(length(allowList), length(spFreq)) = struct('pars', []);
data.thresholds = nan(length(allowList), length(spFreq), length(targetPerf));
data.nearThr = nan(length(allowList), length(spFreq), length(targetPerf));
data.nearThrTrials = cell(length(allowList), length(spFreq), length(targetPerf));
data.thrBnd = nan(length(spFreq), length(targetPerf), 2); % 2 for low and high bounds

% initialize structures for analysis by grouping contrasts
data.contrastGroups = cell(length(allowList), length(spFreq));
data.diffConstantContrast(length(allowList), length(spFreq), nContrastGroups) = ...
    struct('dc', [], 'dsq', [], 'dsqSingleSeg', [],...
    'dcBoot', [], 'useFix', []);
data.saccadeRatesContrast(length(allowList), length(spFreq), nContrastGroups) = ...
    struct('nMS', 0, 'nS', 0, 'total', 0, 'useTrials', []);
data.performanceContrastGroups = nan(length(allowList), length(spFreq), nContrastGroups);
data.performanceContrastGroupsSD = nan(length(allowList), length(spFreq), nContrastGroups);
data.diffConstant(length(allowList), length(spFreq), length(targetPerf)) = ...
    struct('dc', [], 'dsq', [], 'dsqSingleSeg', [],...
    'dcBoot', []);

%% loop through EM filters
fh = 0;
for aa = 1:length(allowList)
    allow = allowList{aa};
    %useTrials excludes no-tracks and blinks and no-response,
    %keeps drift only, only use this spatial frequency
    useTrials = ~counts.notracks(:) & ~counts.blinks(:)...
        & data.Response(:) < 3 ...
        & (data.contrast(:) > 0) & (data.contrast(:) <= 127);
    
    switch allow % filter trials based on eye movements
        case 'all'
            useTrials = data.Response(:) < 3 ...
                & (data.spatialFreq(:) == spFreq(ss)); % no eye movement filtering
        case 'saccades'
            useTrials = useTrials & ...
                (counts.saccades(:) | ...
                counts.microsaccades(:) |...
                counts.driftonly(:));
        case 'microsaccades'
            useTrials = useTrials & ...
                (counts.microsaccades(:) |...
                counts.driftonly(:));
        case 'driftonly'
            useTrials = useTrials & ...
                counts.driftonly(:);
        case 'ms_only'
            useTrials = useTrials & ...
                counts.microsaccades(:) &...
                ~counts.driftonly(:);
        case 'ms_saccades_only'
            useTrials = useTrials & ...
                (counts.saccades(:) | counts.microsaccades(:)) &...
                ~counts.driftonly(:);
    end
    
    fixvalid = useTrials(fixtrialid);
    
    %% psychometric curves for each spatial freq
    for ss = 1:length(spFreq)
        useTrialsSPF = (data.spatialFreq(:) == spFreq(ss)) & useTrials;
        % do my own counts for plotting
        lvl = data.michContrastRounded(useTrialsSPF);
        hits = double(data.Correct(useTrialsSPF));
        ulvl = unique(lvl);
        nlvl = sum(bsxfun(@eq, ulvl(:), lvl), 2);
        ucorr = sum(bsxfun(@times, bsxfun(@eq, ulvl(:), lvl), hits), 2) ./ nlvl;
        
        % do the psychometric fitting for these trials - get contrast
        % thresholds
        [thresh1, par, thresh1B, parB] = psyfit(lvl, hits,...
            'Thresh', targetPerf(1), 'Extra', 'Log', 'Boots', nBoots, 'PlotOff');
        thresh2 = invpsyfun( targetPerf(2), par(1), par(2), ...
            0.5, par(3), false, true);
        thresh2B = nan(size(thresh1B));
        for bb = 1:nBoots
            thresh2B(bb) = invpsyfun( targetPerf(2), parB(1, bb), parB(2, bb), ...
                0.5, parB(3, bb), false, true);
        end
        
        data.thresholds(aa, ss, :) = [thresh1, thresh2]; % mich contrast
        data.psyfit(aa, ss).threshB = [thresh1B(:), thresh2B(:)];
        data.psyfit(aa, ss).pars = par;
        data.psyfit(aa, ss).parsB = parB;
        data.psyfit(aa, ss).lvl = lvl;
        data.psyfit(aa, ss).hits = hits;
        data.psyfit(aa, ss).ulvl = ulvl;
        data.psyfit(aa, ss).nlvl = nlvl;
        data.psyfit(aa, ss).ucorr = ucorr;
        
        % count near threshold trials
        for ti = 1:length(targetPerf)
            % define threshold bounds by performance range
            data.thrBnd(ss, ti, 1) = invpsyfun( max(targetPerf(ti)-perfWidth, .51), par(1), par(2), ...
                0.5, par(3), false, true);
            data.thrBnd(ss, ti, 2) = invpsyfun( min(targetPerf(ti)+perfWidth, .99), par(1), par(2), ...
                0.5, par(3), false, true);
            
            near = (ulvl > data.thrBnd(ss, ti, 1)) & ...
                (ulvl < data.thrBnd(ss, ti, 2));
            data.nearThr(aa, ss, ti) = sum(nlvl(near));
            
            data.nearThrTrials{aa, ss, ti} = useTrialsSPF(:) &...
                (data.michContrastRounded(:) > data.thrBnd(ss, ti, 1)) & ...
                (data.michContrastRounded(:) < data.thrBnd(ss, ti, 2));
            
            assert(data.nearThr(aa, ss, ti) == sum(data.nearThrTrials{aa, ss, ti}));
        end
        
        % now start plotting stuff
        x = logspace(-5, 0, 100);
        p = psyfun(x, par(1), par(2), 0.5, par(3), false, true);
        
        if spFreq(ss) == 6
            xl = [-3, -1]; digs = 3;
        elseif spFreq(ss) == 20
            xl = [-2, 0]; digs = 2;
        end
        xticks = linspace(xl(1), xl(end), 5);
        
        fh = fh + 1;
        figure(fh); clf; hold on;
        for ti = 1:length(targetPerf)
            plot([xl(1), log10(data.thresholds(aa, ss, ti))], targetPerf(ti) * ones(2, 1),...
                'k:', 'linewidth', 2);
            plot(log10(data.thresholds(aa, ss, ti)*ones(2, 1)), [0, targetPerf(ti)],...
                [thrCols(ti), ':'], 'linewidth', 2);
            xarea = log10(x(x > data.thrBnd(ss, ti, 1) & x < data.thrBnd(ss, ti, 2)));
            yarea = p(x > data.thrBnd(ss, ti, 1) & x < data.thrBnd(ss, ti, 2));
            area(xarea(:), yarea(:),...
                'FaceColor', thrCols(ti), 'LineStyle', 'none', 'FaceAlpha', .3);
        end
        
        plot(log10(x), p, [spFCols(ss), '-'], 'linewidth', 2);
        scatter(log10(ulvl), ucorr, 100+nlvl * 30, [spFCols(ss), '.']);
        set(gca, 'XTick', xticks, 'XTicklabel', round(10.^xticks, digs));
        xlabel('Michelson Contrast');
        ylabel('proportion correct');
        title(sprintf('%s - %icpd - %s', subj, spFreq(ss), strrep(allow, '_', ' ')));
        xlim(xl);
        standardPlot;
        
        text(xl(1) + .1 * diff(xl), .95, sprintf('N=%i', sum(useTrialsSPF)),...
            'FontSize', 14);
        for ti = 1:length(targetPerf)
            text(log10(data.thresholds(aa, ss, ti)), .4,...
                sprintf('N_%i=%i', ti, data.nearThr(aa, ss, ti)),...
                'HorizontalAlignment', align{ti},...
                'FontSize', 12);
            text(log10(data.thresholds(aa, ss, ti)), .3,...
                sprintf('thr_%i = %1.3f', ti, data.thresholds(aa, ss, ti)),...
                'HorizontalAlignment', align{ti},...
                'FontSize', 12);
        end
        
        print(fullfile(imgPath, sprintf('%s_psyfun_%icpd_%s.png', subj, spFreq(ss), allow)), '-dpng');
        print(fullfile(imgPath, sprintf('%s_psyfun_%icpd_%s.eps', subj, spFreq(ss), allow)), '-depsc');
    end
    
    %% now analyze drift in near thr trials
    for ss = 1:length(spFreq)
        for ti = 1:length(targetPerf)
            near = data.nearThrTrials{aa, ss, ti}; % accounts for spf and allow
            inear = find(near); % trial indices to use in this analysis
            [~, useTrialst] = intersect(fixtrialid, inear);
            fixtemp = fix(useTrialst);
            
            [~, ~, dc, ~, dsq, dsqSingleSeg] = CalculateDiffusionCoef(fixtemp);
            data.diffConstant(aa, ss, ti).dc = dc;
            data.diffConstant(aa, ss, ti).dsq = dsq;
            data.diffConstant(aa, ss, ti).dsqSingleSeg = dsqSingleSeg;
            
            data.diffConstant(aa, ss, ti).dcBoot = nan(nBoots, 1);
            for bb = 1:nBoots
                rs = randsample(length(fixtemp), length(fixtemp), true);
                [~, ~, dc] = CalculateDiffusionCoef(fixtemp(rs));
                data.diffConstant(aa, ss, ti).dcBoot(bb) = dc;
            end
            
            clear fixtemp;
        end
    end
    
    %% plot diffusion constants by threshold
    fh = fh + 1;
    figure(fh); clf; hold on;
    plot([0, 1], data.fixation.diffConstant.dc * [1, 1], 'k-', 'linewidth', 2);
    plot([0, 1], (data.fixation.diffConstant.dc + std(data.fixation.diffConstant.dcBoot)) * [1, 1], 'k--', 'linewidth', 2);
    plot([0, 1], (data.fixation.diffConstant.dc - std(data.fixation.diffConstant.dcBoot)) * [1, 1], 'k--', 'linewidth', 2);
    hd = nan(length(spFreq), 1);
    for ss = 1:length(spFreq)
        hd(ss) = plot(targetPerf-.02*(ss - 1.5)*2, [data.diffConstant(aa, ss, :).dc], [spFCols(ss), 'o-'], ...
            'markersize', 10, 'linewidth', 2);
        errorbar(targetPerf-.02*(ss - 1.5)*2,...
            [mean(data.diffConstant(aa, ss, 1).dcBoot), mean(data.diffConstant(aa, ss, 2).dcBoot)],...
            [std(data.diffConstant(aa, ss, 1).dcBoot), std(data.diffConstant(aa, ss, 2).dcBoot)],...
            [spFCols(ss), '.'], ...
            'markersize', 15, 'linewidth', 2, 'capsize', 12);
        
    end
    xlim([0.55, 1]);
    legend(hd, cellfun(@(x) sprintf('%icpd', x), (num2cell(spFreq)), 'UniformOutput', false),...
        'Location', 'northeast');
    set(gca, 'XTick', targetPerf, 'XTickLabel', {'low contrast', 'high contrast'});
    ylabel('diffusion constant (arcmin^2 / s)');
    standardPlot;
    title(sprintf('%s - %s', subj, strrep(allow, '_', ' ')));
    
    print(fullfile(imgPath, sprintf('%s_drift_diffCon_thr_%s.png', subj, allow)), '-dpng');
    print(fullfile(imgPath, sprintf('%s_drift_diffCon_thr_%s.eps', subj, allow)), '-depsc');
    
    %% now analyze drift & Sacc grouped by contrast value ( 3 groups )
    
    for ss = 1:length(spFreq)
        contrasts = data.michContrast(useTrials(:) & data.spatialFreq(:) == spFreq(ss));
        clims = quantile(contrasts, linspace(0, 1, nContrastGroups+1));
        data.contrastGroups{aa, ss} = clims;
        for ti = 1:nContrastGroups
            % drift diffusion analysis
            useTrialst = (fixMichContrasts(:) >= clims(ti) & fixMichContrasts(:) < clims(ti+1)) & ...
                fixvalid(:) & (fixspf(:) == spFreq(ss));
            fixtemp = fix(useTrialst);
            
            [~, ~, dc, ~, dsq, dsqSingleSeg] = CalculateDiffusionCoef(fixtemp);
            data.diffConstantContrast(aa, ss, ti).dc = dc;
            data.diffConstantContrast(aa, ss, ti).dsq = dsq;
            data.diffConstantContrast(aa, ss, ti).dsqSingleSeg = dsqSingleSeg;
            data.diffConstantContrast(aa, ss, ti).useFix = useTrialst;
            
            data.diffConstantContrast(aa, ss, ti).dcBoot = nan(nBoots, 1);
            for bb = 1:nBoots
                rs = randsample(length(fixtemp), length(fixtemp), true);
                [~, ~, dc] = CalculateDiffusionCoef(fixtemp(rs));
                data.diffConstantContrast(aa, ss, ti).dcBoot(bb) = dc;
            end
            
            clear fixtemp;
            
            % saccade and microsaccade rates by contrast level
            useTrialst = useTrials(:) & (data.spatialFreq(:) == spFreq(ss)) & ...
                (data.michContrast(:) >= clims(ti) & data.michContrast(:) < clims(ti+1));
            allMSTimes = cat(2, microsaccadeTimes{useTrialst});
            allSTimes = cat(2, saccadeTimes{useTrialst});
            
            nMS = sum(allMSTimes > -diffConDuration & allMSTimes < 0);
            nS = sum(allSTimes > -diffConDuration & allSTimes < 0);
            
            % number of saccades during plat period
            data.saccadeRatesContrast(aa, ss, ti).nMS = nMS / (diffConDuration * sum(useTrialst)) * 1000;
            data.saccadeRatesContrast(aa, ss, ti).nS =  nS / (diffConDuration * sum(useTrialst)) * 1000;
            data.saccadeRatesContrast(aa, ss, ti).total = (nMS + nS) / (diffConDuration * sum(useTrialst)) * 1000;
            data.saccadeRatesContrast(aa, ss, ti).useTrials = useTrialst;
            
            % what's the average performance in these groups?
            correct = data.Correct(useTrialst);
            n = length(correct);
            p = sum(correct) / n;
            data.performanceContrastGroups(aa, ss, ti) = p;
            data.performanceContrastGroupsSD(aa, ss, ti) =...
                sqrt(p * (1-p) / n) * norminv(1-(.05/2));
        end
    end
    
    %% show how trials are grouped on psyfun
    for ss = 1:length(spFreq)
        
        groupcols = 'gmc';
        par = data.psyfit(aa, ss).pars;
        lvl = data.psyfit(aa, ss).lvl;
        hits = data.psyfit(aa, ss).hits;
        ulvl = data.psyfit(aa, ss).ulvl;
        nlvl = data.psyfit(aa, ss).nlvl;
        ucorr = data.psyfit(aa, ss).ucorr;
        
        x = logspace(-5, 0, 100);
        p = psyfun(x, par(1), par(2), 0.5, par(3), false, true);
        
        if spFreq(ss) == 6
            xl = [-4, -1]; digs = 3;
        elseif spFreq(ss) == 20
            xl = [-1.1, 0]; digs = 2;
        end
        xticks = linspace(xl(1), xl(end), 5);
        
        fh = fh + 1;
        figure(fh); clf; hold on;
        for ti = 1:nContrastGroups
            clims = data.contrastGroups{aa, ss}(ti:ti+1);
            xarea = log10(x(x > clims(1) & x < clims(2)));
            yarea = p(x > clims(1) & x < clims(2));
            area(xarea(:), yarea(:), 'FaceColor', groupcols(ti),...
                'LineStyle', 'none', 'FaceAlpha', .3);
        end
        
        plot(log10(x), p, [spFCols(ss), '-'], 'linewidth', 2);
        scatter(log10(ulvl), ucorr, 100+nlvl * 30, [spFCols(ss), '.']);
        set(gca, 'XTick', xticks, 'XTicklabel', round(10.^xticks, digs));
        xlabel('Michelson Contrast');
        ylabel('proportion correct');
        title(sprintf('%s - %icpd - %s', subj, spFreq(ss), strrep(allow, '_', ' ')));
        xlim(xl);
        standardPlot;
        
        text(xl(1) + .1 * diff(xl), .95, sprintf('N=%i', sum(useTrialsSPF)),...
            'FontSize', 14);
        clims = data.contrastGroups{aa, ss};
        clevels = clims(1:end-1) + diff(clims) / 2;
        for ti = 1:nContrastGroups
            text(log10(clevels(ti)), .4,...
                sprintf('N_%i=%i', ti, sum(data.saccadeRatesContrast(aa, ss, ti).useTrials)),...
                'FontSize', 12);
        end
        
        print(fullfile(imgPath, sprintf('%s_psyfun_%icpd_%s_grouped.png', subj, spFreq(ss), allow)), '-dpng');
        print(fullfile(imgPath, sprintf('%s_psyfun_%icpd_%s_grouped.eps', subj, spFreq(ss), allow)), '-depsc');
    end
    
    %% plot diffusion constants by contrast level
    mmm = nan(nContrastGroups, length(spFreq));
    sem = nan(nContrastGroups, length(spFreq));
    for ss = 1:length(spFreq)
        for kk = 1:nContrastGroups
            mmm(kk, ss) = mean(data.diffConstantContrast(aa, ss, kk).dcBoot);
            sem(kk, ss) = std(data.diffConstantContrast(aa, ss, kk).dcBoot);
        end
    end
    
    fh = fh + 1;
    figure(fh); clf; hold on;
    plot([0.0001, 1], data.fixation.diffConstant.dc * [1, 1], 'k-', 'linewidth', 2);
    plot([0.0001, 1], (data.fixation.diffConstant.dc + std(data.fixation.diffConstant.dcBoot)) * [1, 1], 'k--', 'linewidth', 2);
    plot([0.0001, 1], (data.fixation.diffConstant.dc - std(data.fixation.diffConstant.dcBoot)) * [1, 1], 'k--', 'linewidth', 2);
    
    hd = nan(length(spFreq), 1);
    for ss = 1:length(spFreq)
        clevels = data.contrastGroups{aa, ss}(1:end-1) + diff(data.contrastGroups{aa, ss}) / 2;
        hd(ss) = plot(clevels, [data.diffConstantContrast(aa, ss, :).dc], [spFCols(ss), 'o-'], ...
            'markersize', 10, 'linewidth', 2);
        errorbar(clevels, mmm(:, ss), sem(:, ss),...
            [spFCols(ss), '.'], ...
            'markersize', 15, 'linewidth', 2, 'capsize', 12);
        
    end
    xlim(clims([1, end]));
    legend(hd, cellfun(@(x) sprintf('%icpd', x), (num2cell(spFreq)), 'UniformOutput', false),...
        'Location', 'northwest');
    ylabel('diffusion constant (arcmin^2 / s)');
    xlabel('Michelson Contrast');
    standardPlot;
    set(gca, 'XScale', 'log', 'XTick', .1:.2:1); xlim([10^-3.5, 1]);
    title(sprintf('%s - %s', subj, strrep(allow, '_', ' ')));
    
    print(fullfile(imgPath, sprintf('%s_drift_diffCon_contrast_%s.png', subj, allow)), '-dpng');
    print(fullfile(imgPath, sprintf('%s_drift_diffCon_contrast_%s.eps', subj, allow)), '-depsc');
    
    %% plot diffusion constants versus perf grouped by contrast
    fh = fh + 1;
    figure(fh); clf; hold on;
    plot([0, 1], data.fixation.diffConstant.dc * [1, 1], 'k-', 'linewidth', 2);
    plot([0, 1], (data.fixation.diffConstant.dc + std(data.fixation.diffConstant.dcBoot)) * [1, 1], 'k--', 'linewidth', 2);
    plot([0, 1], (data.fixation.diffConstant.dc - std(data.fixation.diffConstant.dcBoot)) * [1, 1], 'k--', 'linewidth', 2);
    
    hd = nan(length(spFreq), 1);
    for ss = 1:length(spFreq)
        clevels = squeeze(data.performanceContrastGroups(aa, ss, :));
        clevelsSD = squeeze(data.performanceContrastGroupsSD(aa, ss, :));
        hd(ss) = plot(clevels, [data.diffConstantContrast(aa, ss, :).dc], [spFCols(ss), 'o-'], ...
            'markersize', 10, 'linewidth', 2);
        errorbar(clevels, mmm(:, ss), ...
            sem(:, ss), sem(:, ss), clevelsSD, clevelsSD,...
            [spFCols(ss), '.'], ...
            'markersize', 15, 'linewidth', 2, 'capsize', 12);
    end
    xlim(clims([1, end]));
    legend(hd, cellfun(@(x) sprintf('%icpd', x), (num2cell(spFreq)), 'UniformOutput', false),...
        'Location', 'northwest');
    ylabel('diffusion constant (arcmin^2 / s)');
    xlabel('proportion correct (by contrast group)');
    standardPlot;
    set(gca,'XTick', .3:.1:1); xlim([.3, 1.1]);
    title(sprintf('%s - %s', subj, strrep(allow, '_', ' ')));
    
    print(fullfile(imgPath, sprintf('%s_drift_diffCon_perf_%s.png', subj, allow)), '-dpng');
    print(fullfile(imgPath, sprintf('%s_drift_diffCon_perf_%s.eps', subj, allow)), '-depsc');
    
    %% plot microsaccades and saccade rates by contrast level
    fh = fh + 1;
    figure(fh); clf; hold on;
    hd = nan(length(spFreq), 1);
    for ss = 1:length(spFreq)
        clevels = data.contrastGroups{aa, ss}(1:end-1) + diff(data.contrastGroups{aa, ss}) / 2;
        hd(ss) = plot(clevels,...
            [data.saccadeRatesContrast(aa, ss, :).total],...
            [spFCols(ss), 'o-'], ...
            'markersize', 10, 'linewidth', 2);
        plot(clevels,...
            [data.saccadeRatesContrast(aa, ss, :).nMS],...
            [spFCols(ss), 'x--'], ...
            'markersize', 10, 'linewidth', 1);
        plot(clevels,...
            [data.saccadeRatesContrast(aa, ss, :).nS],...
            [spFCols(ss), '^:'], ...
            'markersize', 10, 'linewidth', 1);
    end
    set(gca, 'XScale', 'log', 'XTick', .1:.1:1); xlim([10^-3.5, 1]);
    legend(hd, cellfun(@(x) sprintf('%icpd', x), (num2cell(spFreq)), 'UniformOutput', false),...
        'Location', 'northwest');
    ylabel('saccade rate (Hz)');
    xlabel('LOW contrast \leftarrow grating amplitude \rightarrow HIGH contrast');
    standardPlot;
    title(sprintf('%s - %s', subj, strrep(allow, '_', ' ')));
    
    print(fullfile(imgPath, sprintf('%s_saccRate_contrast_%s.png', subj, allow)), '-dpng');
    print(fullfile(imgPath, sprintf('%s_saccRate_contrast_%s.eps', subj, allow)), '-depsc');
end

%% compare to prediction
plotDiffConstantOnPower(subj, data, allowList, spFCols, imgPath);

%% print tables
printTrialCounts(subj, data, counts, spFreq)
printThresholdTables(subj, data);
printDiffConstantTables(subj, data);

%% save data
assert(strcmp(subj, pptrials{1}.Subject_Name(1:4)));
save(fullfile(savePt, sprintf('data.mat')), 'pptrials', 'data', 'counts', 'fix', 'allowList');