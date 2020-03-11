path = 'C:\Users\Ruccilab\Box\Vis_Dynamics\Data';
% subject = 'A013';%'Nikunj'
% subject = 'Nikunj';
subject = 'A092';
pathtodata = fullfile(path,subject);

fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures';
sub_fig_path = fullfile(fig_path,subject,'NeuralModel');
if (~isfolder(sub_fig_path))
    mkdir(sub_fig_path)
end

bin = 0; %set to 1 if binning
read_Data = 0; %if set to 1 all data is read, if 0 previous mat file is loaded
if (read_Data ==1)
    [data, files] = readdata(pathtodata, CalList(), true); %set to true to save in mat file and then use load. 
else
    load(fullfile(pathtodata, 'results.mat'));
end
[ppt] = preprocessing(data); % this is a standard preprocessing call

% ppt = transpose(osRemoverInEM(ppt));
[valid,counter] = countingTrialsNK(ppt);

%%
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
        
        trials(count).x = ppt{i,1}.x;
        trials(count).y = ppt{i,1}.y;
        trials(count).saccades = ppt{i,1}.saccades;
        trials(count).backgroundImage = ppt{i,1}.backgroundImage;
        
        %Correct for erroenous PEST levels
        if trials(count).contrast == 0 && trials(count).present == 1
            trials(count).present = 0;
        end
        
        if trials(count).present == 0
            trials(count).contrast = 0;
        end
        
        if trials(count).contrast > 0.5
            trials(count).contrast = 0.5;
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

%%
e4 = ([trials.ecc] == 4 & [trials.present] == 1); %Trials with stim present at Ecc 4 and 
p1 = [trials.present] == 1;
% inp = trials(e4);
inp = trials(p1);
% tr_num = find(e4==1); %Store actual trial numbers
tr_num = find(p1==1); %Store actual trial numbers
% inp =999; %which trial
% fname = ppt{inp}.backgroundImage;
for i =1:length(inp)
% for i=14:length(inp)
% for i=20
    fname = inp(i).backgroundImage;
    fileID = fopen(fname);
    A = fread(fileID,[1920 1080],'float');
    noiseImg = transpose(A); %the noise image is transposed because of the way they were created in Matlab for EYERIS
    %make grating on top of it
    imsize_x = 1080;
    imsize_y = 1920;
    M_PI = pi;
    width = 1;
    gaussSD = .2;
    noise_amplitude = 0.5;
    eccentricity =  inp(i).ecc;
    SpatialFreq = inp(i).spFreq;
    signal_amplitude = inp(i).contrast;
    phase = 0;
    pixelAngle = inp(i).pixelAngle; % arcmin per pixel, was 0.8 initially
    
    grat = zeros(imsize_x,imsize_y);
    for ii = 1:imsize_x
        for jj = 1:imsize_y
            r = sqrt ((ii - imsize_x/2)^2 + (jj - imsize_y/2)^2) * pixelAngle / 60;
            if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2))
                scale = sin(2 * M_PI * SpatialFreq * r + phase);
            elseif (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD)
                scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                    (cos(M_PI / gaussSD * (r - eccentricity - width / 2)) + 1) / 2;
            elseif (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD)
                scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                    (cos(M_PI / gaussSD * (r - eccentricity + width / 2)) + 1) / 2;
            else
                scale = 0;
            end
            grat(ii, jj) = scale;
        end
    end
    noise1 = (noiseImg - 0.5) * 2;
    noise2 = noise1 * noise_amplitude;
    
    grat1 = signal_amplitude*grat;
    im1 = noise2+grat1;
    im2 = (im1/2)+0.5;
    
    figure(1)
    imshow(im2);
    
    
    
    X = double(inp(i).x.position) + inp(i).xoffset * inp(i).pixelAngle; %degrees
    Y = double(inp(i).y.position) + inp(i).yoffset * inp(i).pixelAngle; %degrees
    
    on = inp(i).saccOn;
    sacc_off = inp(i).saccOff;
    off = inp(i).stimOff;
    cue = inp(i).cueCtr;
    respTime = inp(i).respTime;
    
    [~,emat] = (min(abs(inp(i).saccades.start - on)));
    if isempty(emat)
        continue
    end
    if ~isempty(emat)
        EMAT_saccOn = inp(i).saccades.start(emat);
        EMAT_saccOff = EMAT_saccOn + inp(i).saccades.duration(emat);
    end
    
    x_stim = X;
    y_stim = Y;
    
    x_stim_p = x_stim/pixelAngle;%Pixels
    y_stim_p = y_stim/pixelAngle;%Pixels
    
    x_stim_p_corr = x_stim_p + (imsize_y/2);
    y_stim_p_corr = y_stim_p + (imsize_x/2);
    
   
    [Fc, Fs, R] = MacaqueRetinaM_Spatial_CK(pixelAngle/60,0,11,11,SpatialFreq);
    
    sp_gray = [];
    size = 5;
    for j=1:length(x_stim_p_corr)%looping over time bins
        gray = [];
        gray = [im2(round(y_stim_p_corr(j))-size:1:round(y_stim_p_corr(j))+size,...
            round(x_stim_p_corr(j))-size:1:round(x_stim_p_corr(j))+size)];%get the gray level image of size 30 in noise+grating image around gaze locations
%         figure(3)
%         imshow(gray)
        %     sp_gray=[sp_gray; dot(gray,(Fc - Fs))];%
        %     sp_gray = [sp_gray;sum(sum(gray.*(Fc-Fs)))];
        sp_gray = [sp_gray;sum(sum(gray.*(Fc+Fs)))];
%         sp_gray = [sp_gray;sum(sum(gray.*R))];
    end

    
    [M_Temporal, K, w] = MacaqueRetinaM_Temporal_BK(100);
    res = [];
    res = conv(sp_gray,M_Temporal','same');
    % ... normalize so that the maximum of the response is 1
    res = res/max(abs(res));

    figure(2);
    % plot(EMAT_saccOn:1:round(off),res1);
    % plot(EMAT_saccOn:1:round(off),res);
    
    plot(res,'LineWidth',2);
%     plot(res(cue-50: off+50),'LineWidth',2);
    
    % plot(round(cue):1:round(respTime),res,'LineWidth',2);
    % plot(res(round(cue):round(respTime)),'LineWidth',2);
    % plot(res,'LineWidth',2);
    % plot(EMAT_saccOn+24:1:round(off),res(25:end));
    % plot(EMAT_saccOn:1:round(off),res);
    hold on
    h_emat_sacOn = line([round(EMAT_saccOn) round(EMAT_saccOn)],[round(max(res)) round(min(res))],'LineStyle','--','Color','green','LineWidth', 2);
    h_sacOn = line([round(on) round(on)],[round(max(res)) round(min(res))],'Color','green','LineWidth', 2);
    h_sacOff = line([round(sacc_off) round(sacc_off)],[round(max(res)) round(min(res))],'Color','blue','LineWidth', 2);
    h_emat_sacOff = line([round(EMAT_saccOff) round(EMAT_saccOff)],[round(max(res)) round(min(res))],'LineStyle','--','Color','blue','LineWidth', 2);
    h_stimOff = line([round(off) round(off)],[round(max(res)) round(min(res))],'Color', 'red','LineWidth', 2);
    leg = legend([h_emat_sacOn h_sacOn h_sacOff h_emat_sacOff h_stimOff],'EMAT_saccOn','EYERIS_saccOn','EYERIS_saccOff','EMAT_saccOff','StimOff');
    set(leg,'Interpreter', 'none')
    ylim([-(max(res(round(EMAT_saccOn):round(off)))) (max(res(round(EMAT_saccOn):round(off))))])
%     xlim([round(EMAT_saccOn)-50 round(off)+50])
    set(gca,'box','off','tickdir','out','FontSize',20)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    xlabel('time [ms]')
    ylabel('Response')
    title(['Trial No. ' num2str(tr_num(i)) ', Ecc - ' num2str(eccentricity) ', SpFreq - ' num2str(SpatialFreq) ',Con - ' num2str(signal_amplitude)]);
%     saveas(figure(2),[sub_fig_path '\response_' num2str(tr_num(i))],'epsc')
    hold off

    %EM plot
    figure(3)
    hold on
    hx_1 =  plot(X, 'Color', [0 0 200] / 255, 'HitTest', 'off'); % blue
    hy_1 =  plot(Y, 'Color', [0 180 0] / 255, 'HitTest', 'off'); % green
    
    h_EYERIS_saccOn_1 = line([on on],[-200 500],'Color','green','LineWidth', 2);
    h_EMAT_saccOn_1 = line([EMAT_saccOn EMAT_saccOn],[-200 500],'LineStyle','--','Color','green','LineWidth', 2);
    h_EMAT_saccOff_1 = line([EMAT_saccOff EMAT_saccOff],[-200 500],'LineStyle','--','Color','blue','LineWidth', 2);  
    
    h_stimOff_1 = line([off off],[-200 500],'Color', 'red','LineWidth', 2);
    
    h_EYERIS_saccOff_1 = line([sacc_off sacc_off],[-200 500],'Color','blue','LineWidth', 2);
    
    leg1 = legend([hx_1 hy_1 h_EMAT_saccOn_1 h_EYERIS_saccOn_1 h_EYERIS_saccOff_1 h_EMAT_saccOff_1 h_stimOff_1],...
        'X', 'Y', 'EMAT_saccOn','EYERIS_saccOn','EYERIS_saccOff','EMAT_saccOff', 'StimOff');
    
    set(leg1,'Interpreter', 'none')
    legend boxoff
    grid on
    xlabel('time [ms]')
    ylabel('arcmin')
    set(gca,'box','off','tickdir','out','FontSize',20)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(['Trial No. ',num2str(tr_num(i))]);
    hold off
    
    pause
    close all
    
end

%%
fname = inp(1).backgroundImage;
% fname = strcat(ppt{1}.backgroundImage(1:end-3),'bmp');
% fname = fullfile('C:\Users\Ruccilab\Box\Vis_Dynamics\Matlab_Code',ppt{inp}.backgroundImage); 
fileID = fopen(fname);
% A = fread(fileID,[1080 1920],'float');
A = fread(fileID,[1920 1080],'float');

noiseImg = transpose(A); %the noise image
% imshow(noiseImg);

%make grating on top of it
imsize_x = 1080;
imsize_y = 1920;

M_PI = pi;
width = 1;
gaussSD = .2;
noise_amplitude = 0.5;

eccentricity =  ppt{inp}.eccentricity;
SpatialFreq =  ppt{inp}.spatialFreq;
if ppt{inp}.present
    signal_amplitude = ppt{inp}.contrast;
else
    signal_amplitude = 0;
end
phase = ppt{inp}.phase;
pixelAngle = ppt{inp}.pixelAngle; % arcmin per pixel, was 0.8 initially

grat = zeros(imsize_x,imsize_y);
for ii = 1:imsize_x
    for jj = 1:imsize_y
        r = sqrt ((ii - imsize_x/2)^2 + (jj - imsize_y/2)^2) * pixelAngle / 60;
        if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2))
            scale = sin(2 * M_PI * SpatialFreq * r + phase);
        elseif (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD)
            scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                (cos(M_PI / gaussSD * (r - eccentricity - width / 2)) + 1) / 2;
        elseif (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD)
            scale = sin(2 * M_PI * SpatialFreq * r  + phase) *...
                (cos(M_PI / gaussSD * (r - eccentricity + width / 2)) + 1) / 2;
        else
            scale = 0;
        end
        grat(ii, jj) = scale;
    end
end

noise1 = (noiseImg - 0.5) * 2;
noise2 = noise1 * noise_amplitude;

grat1 = signal_amplitude*grat;
im1 = noise2+grat1;
im2 = (im1/2)+0.5;

figure(2)
imshow(im2);

%%
% X = double(ppt{inp}.x.position) + ppt{inp}.xoffset * ppt{inp}.pixelAngle + eccentricity*60;
% Y = double(ppt{inp}.y.position) + ppt{inp}.yoffset * ppt{inp}.pixelAngle + eccentricity*60;

% vt = osRemoverInEM(ppt{inp});
X = double(ppt{inp}.x.position) + ppt{inp}.xoffset * ppt{inp}.pixelAngle;
Y = double(ppt{inp}.y.position) + ppt{inp}.yoffset * ppt{inp}.pixelAngle;

% X = double(vt.x.position) + ppt{inp}.xoffset * ppt{inp}.pixelAngle;
% Y = double(vt.y.position) + ppt{inp}.yoffset * ppt{inp}.pixelAngle;

% X = double(ppt{inp}.x.position) ;
% Y = double(ppt{inp}.y.position) ;

on = ppt{inp}.saccOn;
sacc_off = ppt{inp}.saccOff;
off = ppt{inp}.stimOff;
cue = ppt{inp}.cueCtr;
respTime = ppt{inp}.responseTime;

[~,emat] = (min(abs(ppt{inp}.saccades.start - on)));
if ~isempty(emat)
    EMAT_saccOn = ppt{inp}.saccades.start(emat);
    EMAT_saccOff = EMAT_saccOn + ppt{inp}.saccades.duration(emat);
end
%Exposure period
% x_stim = X(EMAT_saccOn:off);
% y_stim = Y(EMAT_saccOn:off);

% x_stim = X(cue:respTime);
% y_stim = Y(cue:respTime);
x_stim = X;
y_stim = Y;

% x_stim = X(0:respTime);
% y_stim = Y(0:respTime);

x_stim_p = x_stim/pixelAngle;%Pixels
y_stim_p = y_stim/pixelAngle;

x_stim_p_corr = x_stim_p + (imsize_y/2);
y_stim_p_corr = y_stim_p + (imsize_x/2);

%get gray level values of eye trace during stimulus presentation
gray = [];
% for i=1:length(x_stim_p_corr)
%     gray = [gray im2(round(y_stim_p_corr(i)),round(x_stim_p_corr(i)))];%get the gray level in noise+grating image at gaze locations
% %     gray = [gray im2(round(x_stim_p_corr(i)),round(y_stim_p_corr(i)))];
% end

[Fc, Fs] = MacaqueRetinaM_Spatial_CK(pixelAngle/60,0,101,101,SpatialFreq);
% [Fc, Fs] = MacaqueRetinaM_Spatial_CK(pixelAngle/60,0,SpatialFreq);

sp_gray = [];
size = 50;
for i=1:length(x_stim_p_corr)%looping over time bins
    gray = [];
    gray = [im2(round(y_stim_p_corr(i))-size:1:round(y_stim_p_corr(i))+size,...
        round(x_stim_p_corr(i))-size:1:round(x_stim_p_corr(i))+size)];%get the gray level image of size 30 in noise+grating image around gaze locations
%     imshow(gray)
%     sp_gray=[sp_gray; dot(gray,(Fc - Fs))];% 
%     sp_gray = [sp_gray;sum(sum(gray.*(Fc-Fs)))];
    sp_gray = [sp_gray;sum(sum(gray.*(Fc+Fs)))];
end

% [M_Temporal, K, w] = MacaqueRetinaM_Temporal_BK(length(x_stim_p_corr));
[M_Temporal, K, w] = MacaqueRetinaM_Temporal_BK(100);
res = [];
res = conv(sp_gray,M_Temporal','same');
% res = conv(sp_gray,M_Temporal');
% for j=1:length(sp_gray)
%     res = [res; conv(sp_gray(j,:),M_Temporal')];
% end
% gray2 = round(gray*255);

% f_sq = @(x) (x).^2;
% res1 = f_sq(res);

% f_hwsq = @(x) (x.*(x>=0)).^2; %Half-wave squaring non-linearity f_sq = @(x) (x).^2;
% f_hw = @(x) (x.*(x>=0)); 
% res2 = f_hw(res);
% res = res2;

figure(1);
% plot(EMAT_saccOn:1:round(off),res1);
% plot(EMAT_saccOn:1:round(off),res);
plot(res,'LineWidth',2);
% plot(round(cue):1:round(respTime),res,'LineWidth',2);
% plot(res(round(cue):round(respTime)),'LineWidth',2);
% plot(res,'LineWidth',2);
% plot(EMAT_saccOn+24:1:round(off),res(25:end));
% plot(EMAT_saccOn:1:round(off),res);
hold on
h_emat_sacOn = line([round(EMAT_saccOn) round(EMAT_saccOn)],[round(max(res)) round(min(res))],'LineStyle','--','Color','blue','LineWidth', 2);
h_sacOn = line([round(on) round(on)],[round(max(res)) round(min(res))],'LineStyle','--','Color','green','LineWidth', 2);
h_sacOff = line([round(sacc_off) round(sacc_off)],[round(max(res)) round(min(res))],'LineStyle','--','Color','k','LineWidth', 2);
h_emat_sacOff = line([round(EMAT_saccOff) round(EMAT_saccOff)],[round(max(res)) round(min(res))],'LineStyle','--','Color','yellow','LineWidth', 2);
h_stimOff = line([round(off) round(off)],[round(max(res)) round(min(res))],'LineStyle','--','Color','red','LineWidth', 2);
leg = legend([h_emat_sacOn h_sacOn h_sacOff h_emat_sacOff h_stimOff],'EMAT_saccOn','EYERIS_saccOn','EYERIS_saccOff','EMAT_saccOff','StimOff');
set(leg,'Interpreter', 'none')
% ylim([0 5000])
% xlim([round(EMAT_saccOn)-10 round(off)+10])
set(gca,'box','off','tickdir','out','FontSize',20)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xlabel('time [ms]')
ylabel('Response')
title(['Trial No. ' num2str(inp) ', Ecc - ' num2str(eccentricity) ', SpFreq - ' num2str(SpatialFreq)]);
saveas(figure(1),[sub_fig_path '\response_' num2str(inp)],'epsc')


%%
% figure(1);
% plot(gray2)
% set(gca,'box','off','tickdir','out','FontSize',30)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% xlabel('Time (ms)')
% ylabel('Gray level')
% title('Retinal input')
% saveas(figure(1),[sub_fig_path '\retinalinput'],'epsc');
% 
% %Single neuron
% spatial_FR = gray2.(Fc - Fs);
% w= conv(gray2,MacaqueRetinaM_Temporal_BK);
% 
% figure(2);
% plot(w,'linewidth', 2)
% set(gca,'box','off','tickdir','out','FontSize',30)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% xlabel('Time (ms)')
% ylabel('Response')
% title('Neuronal response M-cell')
% saveas(figure(2),[sub_fig_path '\responseM'],'epsc');
% 
% w2= conv(gray2,MacaqueRetinaP_Temporal_BK);
% 
% figure(3);
% plot(w2,'linewidth', 2)
% set(gca,'box','off','tickdir','out','FontSize',30)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% xlabel('Time (ms)')
% ylabel('Response')
% title('Neuronal response P-cell')
% saveas(figure(3),[sub_fig_path '\responseP'],'epsc');