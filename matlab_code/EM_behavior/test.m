%% distribution of saccade directions (polar)
fh = fh + 1;
figure(fh); clf;
h = nan(length(plottasks), 2);
for ti = plottasks
    tmpc = out.(sprintf('dir_%s_hist_bins', tasks{ti}));
    tmph = out.(sprintf('dir_%s_hist', tasks{ti}));
    for EYE = 1:2
        c = tmpc{EYE};
        n1 = tmph{EYE};
        h(ti, EYE) = polar([c, c(1)], [n1, n1(1)]);
        set(h(ti, EYE), 'Color', cols4{ti, EYE});
        hold on;
    end
end
set(h, 'linewidth', 2);
xlabel('sacc direction prob');
legend(h(:, 1), tasks);
mysavefig('ms-dir-distribution');

%% nested function to plot 1d distributions (not polar) across tasks & eye
    function plot_1d_dist(prop)
        hh = nan(length(plottasks), 2);
        for tii = plottasks
            tmpcc = out.(sprintf('%s_%s_hist_bins', prop, tasks{tii}));
            tmphh = out.(sprintf('%s_%s_hist', prop, tasks{tii}));
            for EYEi = 1:2
                if length(tmphh) < EYEi
                    break; 
                end
                hh(tii, EYEi) = plot(tmpcc{EYEi}(1:end-1), tmphh{EYEi}(1:end-1),...
                    'linewidth', 2, 'DisplayName', tasks{tii},...
                    'Color', cols4{tii, EYEi});
            end
        end
        
        for tii = plottasks
            tmpms = out.(sprintf('%s_%s_MS', prop, tasks{tii}));
            for EYEi = 1:2
                if length(tmpms) < EYEi
                    break; 
                end
                vertLineThrough(tmpms{EYEi}(1), cols4{tii, EYEi}, gca, '--');
            end
        end
        
        standardPlot;
        legend(hh(:, 1), tasks(plottasks), 'Location', 'northeast', 'FontSize', LFONTSIZE);
    end