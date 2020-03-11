fig_path='C:\Users\Ruccilab\Box\Vis_Dynamics\Figures';
load_BOOT_mat =1;%set to 1 if loading from saved matfile
subjects = ["Nikunj","A036","A013","A092",];
% subjects = ["Nikunj"];
for si = 1:length(subjects)

    sub_fig_path = fullfile(fig_path,subjects(si),'psychfit');
    % sub_fig_path = fullfile(fig_path,subject,'land30');
    if (~isfolder(sub_fig_path))
        mkdir(sub_fig_path)
    end
    if load_BOOT_mat ==1
        matfilename = strcat(subjects(si),'_plot_BOOT.mat');
        load(matfilename,'plotting')
    else
        matfilename = strcat(subjects(si),'_plot.mat');
        load(matfilename,'plotting')
    end
    ecc_levels = unique([plotting.ecc]);
    contrast = [plotting.contrast];
    pres_Time = [plotting.presTime];
    sp_Freq = [plotting.spFreq];
    pres_levels= unique([plotting.presTime]);
    freq_levels= unique([plotting.spFreq]);
    
    
    
    
    %Sensitivity/Threshold summary plots
    %Pres time vs sensitivity (2 cpd)
    fh=1;
    figure(fh)
    alpha = [];
    % for c=1:length(freq_levels)
    %     for p = 1:length(pres_levels)
    for i=1:length(ecc_levels)
        In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
        In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
        t1 = plotting(In1).thresh;
        t2 = plotting(In2).thresh;
%         outliers = find(plotting(In1).threshB>1);
        plotting(In1).threshB(plotting(In1).threshB>1) = [];
        plotting(In1).threshB(plotting(In1).threshB<0) = [];
        std1 = std(plotting(In1).threshB);
        std2 = std(plotting(In2).threshB);
        
        stddev = [1/std1 1/std2];
        if t1>0.6
            t1 = 0.6;
        elseif t1<0
            t1 = 0.01;
        end
        if t2>0.6
            t2 = 0.6;
        elseif t2<0
            t2 = 0.01;
        end
        sens = [1/t1 1/t2];
        alpha = [alpha (sens(2) - sens(1))/(450)];
        formatSpec = "(sens/ms)";
        %             if sens>0.5
        %                 sens = 0.5;
        %             elseif sens<0
        %                 sens = 0;
        %             end
        if ecc_levels(i)==0
            plot([50 500],sens,'-s','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            %                 errorbar([50 500],sens,stddev,'-s','linewidth', 4,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            txt = ['\angle' alpha(i) formatSpec];
            %                 txt = ['\angle',alpha(i)];
            text(50,sens(1),txt,'FontSize',20)
            hold on
        elseif ecc_levels(i)==4
            plot([50 500],sens,'-v','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            %                 errorbar([50 500],sens,stddev,'-v','linewidth', 4,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            txt = ['\angle' alpha(i) formatSpec];
            %                 txt = ['\angle',alpha(i)];
            text(50,sens(1),txt,'FontSize',20)
            hold on
        elseif ecc_levels(i)==8
            plot([50 500],sens,'-o','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            %                 errorbar([50 500],sens,stddev,'-o','linewidth', 4,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            txt = ['\angle' alpha(i) formatSpec];
            %                 txt = ['\angle',alpha(i)];
            text(50,sens(1),txt,'FontSize',20)
            hold on
        end
    end
    
    hold off
    set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    legend
    xticks([50 500])
    xlim([-50 600])
    xlabel('Presentation time')
    ylabel('Sensitivity')
    title('Low SpFreq')
%     saveas(figure(fh),[sub_fig_path '\Sens_Pres_2'],'epsc');
    saveas(figure(fh),fullfile(sub_fig_path,'\Sens_Pres_2'),'epsc');
    
    %Pres time vs sensitivity (10 cpd)
    
    alpha = [];
    fh=fh+1;
    figure(fh)
    % y_tick=[];
    % for c=1:length(freq_levels)
    %     for p = 1:length(pres_levels)
    for i=1:length(ecc_levels)
        In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
        In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
        t1 = plotting(In1).thresh;
        t2 = plotting(In2).thresh;
        plotting(In1).threshB(plotting(In1).threshB>1) = [];
        plotting(In1).threshB(plotting(In1).threshB<0) = [];
        std1 = std(plotting(In1).threshB);
        std2 = std(plotting(In2).threshB);
        stddev = [1/std1 1/std2];
        if t1>0.6
            t1 = 0.6;
        elseif t1<0
            t1 = 0.01;
        end
        if t2>0.6
            t2 = 0.6;
        elseif t2<0
            t2 = 0.01;
        end
        sens = [1/t1 1/t2];
        %             alpha = [alpha atand((sens(2) - sens(1))/(450))];
        alpha = [alpha (sens(2) - sens(1))/(450)];
        formatSpec = "(sens/ms)";
        %             if ecc_levels(i)==0 && freq_levels(c)==2
        %                 plot([50 500],th,'-s','linewidth', 2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
        %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
        %                 hold on
        %             elseif ecc_levels(i)==4 && freq_levels(c)==2
        %                 plot([50 500],th,'-v','linewidth', 2,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
        %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
        %                 hold on
        %             elseif ecc_levels(i)==8 && freq_levels(c)==2
        %                 plot([50 500],th,'-o','linewidth', 2,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
        %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
        %                 hold on
        if ecc_levels(i)==0
            plot([50 500],sens,'-s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            %                 errorbar([50 500],sens,stddev,'-s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            %                 txt = sprintf(formatSpec,alpha(i));
            txt = ['\angle' alpha(i) formatSpec];
            %                 txt = ['\angle',alpha(i)];
            text(50,sens(1),txt,'FontSize',20)
            hold on
        elseif ecc_levels(i)==4
            %                 errorbar([50 500],sens,stddev,'-v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            plot([50 500],sens,'-v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            txt = ['\angle' alpha(i) formatSpec];
            text(50,sens(1),txt,'FontSize',20)
            hold on
        elseif ecc_levels(i)==8
            %                 errorbar([50 500],sens,stddev,'-o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            plot([50 500],sens,'-o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            txt = ['\angle' alpha(i) formatSpec];
            text(50,sens(1),txt,'FontSize',20)
            hold on
        end
        %             y_tick = [y_tick th];
    end
    %     end
    %
    % end
    hold off
    set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    legend
    xticks([50 500])
    xlim([-50 600])
%     ylim([0 0.6])
    % yticks(sort(y_tick))
    xlabel('Presentation time')
    ylabel('Sensitivity')
    title('High SpFreq')
%     saveas(figure(fh),[sub_fig_path '\Sens_Pres_10'],'epsc');
    saveas(figure(fh),fullfile(sub_fig_path,'\Sens_Pres_10'),'epsc');
    
    %Pres time vs threshold (2 cpd)
    fh=fh+1;
    figure(fh)
    % for c=1:length(freq_levels)
    %     for p = 1:length(pres_levels)
    for i=1:length(ecc_levels)
        In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
        In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==2); %indices
        t1 = plotting(In1).thresh;
        t2 = plotting(In2).thresh;
        plotting(In1).threshB(plotting(In1).threshB>1) = [];
        plotting(In1).threshB(plotting(In1).threshB<0) = [];
        plotting(In2).threshB(plotting(In2).threshB>1) = [];
        plotting(In2).threshB(plotting(In2).threshB<0) = [];
        std1 = std(plotting(In1).threshB);
        std2 = std(plotting(In2).threshB);
        stddev = [std1 std2];
        if t1>0.6
            t1 = 0.6;
        elseif t1<0
            t1 = 0.01;
        end
        if t2>0.6
            t2 = 0.6;
        elseif t2<0
            t2 = 0.01;
        end
        thresh = [t1 t2];
        if ecc_levels(i)==0
            errorbar([50 500],thresh,stddev,'--s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            hold on
        elseif ecc_levels(i)==4
            errorbar([50 500],thresh,stddev,'--v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            hold on
        elseif ecc_levels(i)==8
            errorbar([50 500],thresh,stddev,'--o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 2'])
            hold on
            %             elseif ecc_levels(i)==0 && freq_levels(c)==10
            %                 plot([50 500],th,'--s','linewidth', 2,'Color','k','MarkerEdgeColor','k',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
            %                 hold on
            %             elseif ecc_levels(i)==4 && freq_levels(c)==10
            %                 plot([50 500],th,'--v','linewidth', 2,'Color','r','MarkerEdgeColor','r',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
            %                 hold on
            %             elseif ecc_levels(i)==8 && freq_levels(c)==10
            %                 plot([50 500],th,'--o','linewidth', 2,'Color','b','MarkerEdgeColor','b',...
            %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
            %                 hold on
        end
    end
    %     end
    
    % end
    hold off
    set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    legend
    xticks([50 500])
    xlim([-50 600])
    ylim([0 0.6])
    xlabel('Presentation time')
    ylabel('Threshold')
    title('Low SpFreq')
%     saveas(figure(fh),[sub_fig_path '\Thresh_Pres_2'],'epsc');
    saveas(figure(fh),fullfile(sub_fig_path,'\Thresh_Pres_2'),'epsc');
    
    %Pres time vs threshold (10 cpd)
    fh=fh+1;
    figure(fh)
    % y_tick=[];
    % for c=1:length(freq_levels)
    %     for p = 1:length(pres_levels)
    for i=1:length(ecc_levels)
        In1=([plotting.presTime]==50 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
        In2=([plotting.presTime]==500 & [plotting.ecc]==ecc_levels(i) & [plotting.spFreq]==10); %indices
        t1 = plotting(In1).thresh;
        t2 = plotting(In2).thresh;
        plotting(In1).threshB(plotting(In1).threshB>1) = [];
        plotting(In1).threshB(plotting(In1).threshB<0) = [];
        plotting(In2).threshB(plotting(In2).threshB>1) = [];
        plotting(In2).threshB(plotting(In2).threshB<0) = [];
        std1 = std(plotting(In1).threshB);
        std2 = std(plotting(In2).threshB);
        stddev = [std1 std2];
        if t1>0.6
            t1 = 0.6;
        elseif t1<0
            t1 = 0.01;
        end
        if t2>0.6
            t2 = 0.6;
        elseif t2<0
            t2 = 0.01;
        end
        thresh = [t1 t2];
        
        %             if ecc_levels(i)==0 && freq_levels(c)==2
        %                 plot([50 500],th,'-s','linewidth', 2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k',...
        %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
        %                 hold on
        %             elseif ecc_levels(i)==4 && freq_levels(c)==2
        %                 plot([50 500],th,'-v','linewidth', 2,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r',...
        %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
        %                 hold on
        %             elseif ecc_levels(i)==8 && freq_levels(c)==2
        %                 plot([50 500],th,'-o','linewidth', 2,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b',...
        %                     'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - ',num2str(freq_levels(c))])
        %                 hold on
        if ecc_levels(i)==0
            errorbar([50 500],thresh,stddev,'--s','linewidth', 4,'Color','k','MarkerEdgeColor','k',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            hold on
        elseif ecc_levels(i)==4
            errorbar([50 500],thresh,stddev,'--v','linewidth', 4,'Color','r','MarkerEdgeColor','r',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            hold on
        elseif ecc_levels(i)==8
            errorbar([50 500],thresh,stddev,'--o','linewidth', 4,'Color','b','MarkerEdgeColor','b',...
                'DisplayName',['Ecc - ',num2str(ecc_levels(i)),',SpFreq - 10'])
            hold on
        end
        %             y_tick = [y_tick th];
    end
    %     end
    %
    % end
    hold off
    set(gca,'box','off','tickdir','out','FontSize',30,'YScale', 'log')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    legend
    xticks([50 500])
    xlim([-50 600])
    ylim([0 0.6])
    % yticks(sort(y_tick))
    xlabel('Presentation time')
    ylabel('Threshold')
    title('High SpFreq')
%     saveas(figure(fh),[sub_fig_path '\Thresh_Pres_10'],'epsc');
    saveas(figure(fh),fullfile(sub_fig_path,'\Thresh_Pres_10'),'epsc');
    %% Radial Colorbar plot of thresholds
    
%     imsize_x = 1080;
%     imsize_y = 1920;
%     pixelAngle = 1.06; % arcmin per pixel, was 0.8 initially
%     width = 1;
%     fh=fh+1;
%     fig_count =fh;
%     for p = 1:length(pres_levels) %different presentation times
%         for c=1:length(freq_levels)
%             grat = zeros(imsize_x,imsize_y);
%             for ii = 1:imsize_x
%                 for jj = 1:imsize_y
%                     r = sqrt ((ii - imsize_x/2)^2 + (jj - imsize_y/2)^2)* pixelAngle / 60;
%                     if ((r >= ecc_levels(1) - width/2) && (r <= ecc_levels(1) + width/2))
%                         In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(1) & [plotting.spFreq]==freq_levels(c)); %indices
%                         scale = plotting(In1).thresh;
%                     elseif ((r >= ecc_levels(2) - width/2) && (r <= ecc_levels(2) + width/2))
%                         In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(2) & [plotting.spFreq]==freq_levels(c)); %indices
%                         scale = plotting(In1).thresh;
%                     elseif ((r >= ecc_levels(3) - width/2) && (r <= ecc_levels(3) + width/2))
%                         In1=([plotting.presTime]==pres_levels(p) & [plotting.ecc]==ecc_levels(3) & [plotting.spFreq]==freq_levels(c)); %indices
%                         scale = plotting(In1).thresh;
%                     else
%                         scale = 0;
%                     end
%                     grat(ii, jj) = scale;
%                 end
%             end
%             figure(fig_count)
%             imagesc(grat);
%             set (gca, 'CLim',[0 0.5])
%             colorbar
%             title(['Pres. - ',num2str(pres_levels(p)),' ,Freq. - ',num2str(freq_levels(c))])
% %             saveas(figure(fig_count),[sub_fig_path '\Color_' ['Pres_',num2str(pres_levels(p)),'_Freq_',num2str(freq_levels(c))]],'epsc');
%             saveas(figure(fig_count),fullfile(sub_fig_path,['Color_Pres_',num2str(pres_levels(p)),'_Freq_',num2str(freq_levels(c))]),'epsc');
%             fig_count =fig_count+1;
%         end
%     end
    
   close all
end

%%
%Save threshold summary table

% thresh_T = table;
% formatSpec = 'Ecc-%d, Pres-%d, SpFreq-%d';
% dummy_cell = {};
% for i =1:12
%     name = sprintf(formatSpec,plotting(i).ecc,plotting(i).presTime,plotting(i).spFreq);
%     row_1 = {name,num2str(round(plotting(i).thresh,2))};
%     dummy_cell{i,1} = row_1;
%     colnames = {'Condition';'Threshold'};
%     dummy = cell2table(dummy_cell{i,1},'VariableNames',colnames);
%     thresh_T = [thresh_T; dummy];    
% end
% writetable(thresh_T,[sub_fig_path '\psychfit.csv'],'WriteVariableNames',true);

