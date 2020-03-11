function displayTrialSequence(vt, trials)

BackgroundColor = [0.85 0.85 0.85];
MainFigure = figure('Visible', 'on', ...
    'Name', 'Display Trial', ...
    'Position', [400 100 1200 830], ...
    'NumberTitle', 'off', ...
    'Toolbar', 'figure', ...
    'Resize', 'off', 'Color', BackgroundColor);
set(MainFigure, 'PaperPositionMode','auto');
set(MainFigure, 'InvertHardcopy','off');
set(MainFigure, 'Units','pixels');
[pxAngle] = CalculatePxAngle(vt{1}.Xres, 1550, 600);
vt = osRemoverInEM(vt);

DelayAfterTOn = 450;


for idx = 1:length(trials)
    ii = trials(idx);
    X = double(vt{ii}.x.position) + vt{ii}.xoffset * pxAngle;
    Y = double(vt{ii}.y.position) + vt{ii}.yoffset * pxAngle;
    
    FixationON = vt{ii}.TimeFixationON;
    FixationOFF = FixationON + vt{ii}.FixationTime; 
    SaccCueON = FixationOFF + round(vt{ii}.BeepDelayTime);
    SaccCueOFF = SaccCueON + vt{ii}.SaccCueTime;
    TargetON = vt{ii}.TimeTargetON;    
    TargetOFF = TargetON + vt{ii}.TargetTime;
    RespCueON = vt{ii}.TimeCueON;
    RespCueOFF = RespCueON + vt{ii}.CueTime; 
    ResponseTime = vt{ii}.ResponseTime;
       
    hold on
    Limits = [-50,50];
    
    % cue
    h_c = fill([RespCueON RespCueON ...
        RespCueON + RespCueOFF-RespCueON ...
        RespCueON + RespCueOFF-RespCueON], ...
        [Limits(1) (Limits(2) - 0.5) (Limits(2) - 0.5) Limits(1)], ...
        [255 175 100]/255, ...
        'EdgeColor', [255 175 100]/255, ...
        'LineWidth', 2, 'Clipping', 'On');
    
    % target
    h_t = fill([TargetON TargetON ...
        TargetON + TargetOFF-TargetON ...
        TargetON + TargetOFF-TargetON], ...
        [Limits(1) (Limits(2) - 0.5) (Limits(2) - 0.5) Limits(1)], ...
        [175 230 100]/255, ...
        'EdgeColor', [175 230 100]/255, ...
        'LineWidth', 2, 'Clipping', 'On');
    
    % saccade cue
    h_s = fill([SaccCueON SaccCueON ...
        SaccCueON + SaccCueOFF-SaccCueON ...
        SaccCueON + SaccCueOFF-SaccCueON], ...
        [Limits(1) (Limits(2) - 0.5) (Limits(2) - 0.5) Limits(1)], ...
        [175 80 20]/255, ...
        'EdgeColor', [175 80 20]/255, ...
        'LineWidth', 2, 'Clipping', 'On');
    
    % response
    h_r = fill([ResponseTime ResponseTime ...
        ResponseTime + ResponseTime+50-ResponseTime ...
        ResponseTime + ResponseTime+50-ResponseTime], ...
        [Limits(1) (Limits(2) - 0.5) (Limits(2) - 0.5) Limits(1)], ...
        [70 100 20]/255, ...
        'EdgeColor', [70 100 20]/255, ...
        'LineWidth', 1, 'Clipping', 'On');
    
        % valid ms time interval
    h_validMs = fill([TargetON + 80, TargetON + 80, ...
        RespCueON, ...
        RespCueON], ...
        [Limits(1), (Limits(2) - 0.5), (Limits(2) - 0.5), Limits(1)], ...
        'b', 'LineWidth', 2, 'EdgeColor', 'b', ...
        'Clipping', 'On', 'FaceAlpha', .25);
    
%     h_noMs = fill([TargetON - 100, TargetON - 100, ...
%         TargetON + DelayAfterTOn, ...
%         TargetON + DelayAfterTOn], ...
%         [Limits(1), (Limits(2) - 0.5), (Limits(2) - 0.5), Limits(1)], ...
%         'm', 'LineWidth', 2, 'EdgeColor', 'm', ...
%         'Clipping', 'On', 'FaceAlpha', .25);
    
    h_clearMs = fill([SaccCueON - 80, SaccCueON - 80, ...
        TargetON, TargetON], ...
        [Limits(1), (Limits(2) - 0.5), (Limits(2) - 0.5), Limits(1)], ...
        'r', 'LineWidth', 2, 'EdgeColor', 'r', ...
        'Clipping', 'On', 'FaceAlpha', .25);
    
    hx =  plot(X, 'Color', [0 0 200] / 255, 'HitTest', 'off'); % blue
    hy =  plot(Y, 'Color', [0 180 0] / 255, 'HitTest', 'off'); % green
    ylim([-35,35])
%     if vt{ii}.SaccCueType == 1
%         TposR = plot([0 5000], [vt{ii}.TargetOffsetpx * pxAngle, vt{ii}.TargetOffsetpx * pxAngle], 'g', 'LineWidth', 3);
%         TposL = plot([0 5000], [-vt{ii}.TargetOffsetpx * pxAngle, -vt{ii}.TargetOffsetpx * pxAngle], 'k', 'LineWidth', 1);
%     elseif vt{ii}.SaccCueType == 3
%         TposR = plot([0 5000], [vt{ii}.TargetOffsetpx * pxAngle, vt{ii}.TargetOffsetpx * pxAngle], 'b', 'LineWidth', 3);
%         TposL = plot([0 5000], [-vt{ii}.TargetOffsetpx * pxAngle, -vt{ii}.TargetOffsetpx * pxAngle], 'k', 'LineWidth', 1);
%     elseif vt{ii}.SaccCueType == 5
%         TposR = plot([0 5000], [vt{ii}.TargetOffsetpx * pxAngle, vt{ii}.TargetOffsetpx * pxAngle], 'k', 'LineWidth', 1);
%         TposL = plot([0 5000], [-vt{ii}.TargetOffsetpx * pxAngle, -vt{ii}.TargetOffsetpx * pxAngle], 'g', 'LineWidth', 3);
%     elseif vt{ii}.SaccCueType == 7
%         TposR = plot([0 5000], [vt{ii}.TargetOffsetpx * pxAngle, vt{ii}.TargetOffsetpx * pxAngle], 'k', 'LineWidth', 1);
%         TposL = plot([0 5000], [-vt{ii}.TargetOffsetpx * pxAngle, -vt{ii}.TargetOffsetpx * pxAngle], 'b', 'LineWidth', 3);
%     end
    
%     % sacc cue direction
%     if vt{ii}.SaccCueType == 5
%         plot(1,vt{ii}.TargetOffsetpx*pxAngle, '*', 'LineWidth', 10)
%     elseif vt{ii}.SaccCueType == 2
%         plot(1,-vt{ii}.TargetOffsetpx*pxAngle, '*', 'LineWidth', 10)
%     elseif vt{ii}.SaccCueType == 0
%         plot(1,vt{ii}.TargetOffsetpx*pxAngle, '*', 'LineWidth', 10)
%         plot(1,-vt{ii}.TargetOffsetpx*pxAngle, '*', 'LineWidth', 10)
%     end
%     

    legend boxoff
    legend([hx hy h_t h_s h_c h_r h_validMs h_clearMs], 'X', 'Y', 'Target','Sacc. Cue','Resp. Cue', 'Response', ...
        'Period of Interest','No Ms Allowed', 'Location','EastOutside','Orientation','Vertical')
    legend boxoff
    set(gca,'FontSize', 15, 'YTick', -20:20:20, 'XTick', 0:500:7000)
    grid on
    xlabel('time [ms]')
    ylabel('arcmin')
    
    %xlim([1500 2500])
    %     if (vt{ii}.TrialType==1 & sign(vt{ii}.X_TargetOffset)==1)
    %     title(sprintf('Valid trial, cue to the right, RT = %.2f', vt{ii}.RT))
    %     elseif  (vt{ii}.TrialType==1 & (sign(vt{ii}.X_TargetOffset)==-1))
    %      title(sprintf('Valid trial, cue to the left, RT = %.2f', vt{ii}.RT))
    %     elseif  (vt{ii}.TrialType==0 & (sign(vt{ii}.X_TargetOffset)==-1))
    %      title(sprintf('Invalid trial, cue to the left, RT = %.2f', vt{ii}.RT))
    %     elseif (vt{ii}.TrialType==0 & sign(vt{ii}.X_TargetOffset)==1)
    %     title(sprintf('Invalid trial, cue to the right, RT = %.2f', vt{ii}.RT))
    %     elseif  (vt{ii}.TrialType==2 & (sign(vt{ii}.X_TargetOffset)==-1))
    %         title(sprintf('Neutral trial, target to the left, RT = %.2f', vt{ii}.RT))
    %     elseif (vt{ii}.TrialType==2 & sign(vt{ii}.X_TargetOffset)==1)
    %         title(sprintf('Neutral trial, target to the right, RT = %.2f', vt{ii}.RT))
    %     elseif  (vt{ii}.TrialType==3 & (sign(vt{ii}.X_TargetOffset)==-1))
    %         title(sprintf('Catch trial, cue to the left, RT = %.2f', vt{ii}.RT))
    %     elseif (vt{ii}.TrialType==3 & sign(vt{ii}.X_TargetOffset)==1)
    %         title(sprintf('Catch trial, cue to the right, RT = %.2f', vt{ii}.RT))
    %     end
    
    input ''
    
end


end