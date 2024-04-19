%% 
% Load data

clear variables

monkey = 'P';
load(['\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results_StatsFinal'])


All_day_index = cellfun(@(x) ones(size(x,3),1), signif_units, 'UniformOutput', false);
for d = 1:length(signif_units)
    All_day_index{d} = d*ones(size(signif_units{d},3),1);
end
All_day_index = cat(1,All_day_index{:});
    
    
All_avg_speed_init = avg_speed_init(All_day_index);
All_avg_speed_cor  = avg_speed_cor(All_day_index);

All_signif_units = squeeze(cat(3,signif_units{:}));
All_Coefficients_init = cat(2,best_Coefficients_init{:});
All_Coefficients_cor = cat(2,best_Coefficients_cor{:});

All_Coefficients_init1 = squeeze(cat(3,best_Coefficients_init1{:}));
All_Coefficients_init2 = squeeze(cat(3,best_Coefficients_init2{:}));
All_Coefficients_cor1 = squeeze(cat(3,best_Coefficients_cor1{:}));
All_Coefficients_cor2 = squeeze(cat(3,best_Coefficients_cor2{:}));


%Velocity regression coefficients depth of modulation

vel_Modulation_init = 2*sqrt(sum(All_Coefficients_init(4:5,:).^2,1));
vel_Modulation_cor  = 2*sqrt(sum(All_Coefficients_cor(4:5,:).^2,1));

vel_Modulation_init1 = squeeze( 2*sqrt(sum(All_Coefficients_init1(4:5,:,:).^2,1)) );
vel_Modulation_init2 = squeeze( 2*sqrt(sum(All_Coefficients_init2(4:5,:,:).^2,1)) );
vel_Modulation_cor1  = squeeze( 2*sqrt(sum(All_Coefficients_cor1(4:5,:,:).^2,1)) );
vel_Modulation_cor2  = squeeze( 2*sqrt(sum(All_Coefficients_cor2(4:5,:,:).^2,1)) );


    
%Data regression coefficients depth of modulation
data_Modulation_init = 2*All_avg_speed_init.*vel_Modulation_init;
data_Modulation_cor  = 2*All_avg_speed_cor.*vel_Modulation_cor;

%%
% Graphs

init_color = [0,0,0.7];
cor_color = [0.7,0,0];


% Calculate vel_mod correlation coefficient
[vel_mod_correlation_matrix,vmP,vmRL,vmRU] = corrcoef(vel_Modulation_init(All_signif_units), vel_Modulation_cor(All_signif_units));
vel_mod_correlation = vel_mod_correlation_matrix(1, 2);
disp(['Correlation init vs. corr speed modifed: ' num2str(vel_mod_correlation), ' 95% CI: ' num2str(vmRL(1,2)) '-' num2str(vmRU(1,2))])
%% Correlation init vs. corr speed modifed: 0.47165 95% CI: 0.4136-0.52588

figure
scatter(vel_Modulation_init(All_signif_units), vel_Modulation_cor(All_signif_units), 'MarkerEdgeColor', 'k')
hold on
line([0 45], [0 90], 'LineStyle', '--', 'Color', 'k')
line([0 90], [0 90], 'LineStyle', '--', 'Color', 'k')
line([0 90], [0 45], 'LineStyle', '--', 'Color', 'k')
xlim([-0.1,90])
ylim([-0.1,90])
text(75,32,'0.5x', 'FontSize', 12)
text(70,65,'1x', 'FontSize', 12)
text(39,75,'2x', 'FontSize', 12)
axis square
xlabel('Init. Depth of Mod. (Hz/Vel)', 'Color', init_color)
ylabel('Corr. Depth of Mod. (Hz/Vel)', 'Color', cor_color)
set(gca, 'FontSize', 14)
% set(gcf,'Position',[0 0 600 400]);
% title(['Monkey ' monkey])
% print(gcf, ['./Schwartze_Paper1/figure_graphics/' monkey '_figure5_Vel_DOM_p'], '-dpng')
% if monkey == 'P'
%     fig2svg('./figure_graphics/figure5d.svg', gcf)
% elseif monkey == 'Q'
%     fig2svg('./figure_graphics/figure5e.svg', gcf)
% end

figure
scatter(data_Modulation_init(All_signif_units), data_Modulation_cor(All_signif_units), 'MarkerEdgeColor', 'k')
hold on
line([0 150], [0 160], 'LineStyle', '--', 'Color', 'k')
xlim([-0.1,150])
ylim([-0.1,150])
axis square
xlabel('Init. FR Mod. (Hz)', 'Color', init_color)
ylabel('Corr. FR Mod. (Hz)', 'Color', cor_color)
set(gca, 'FontSize', 14)
% set(gcf,'Position',[0 0 600 400]);
% title(['Monkey ' monkey])
% print(gcf, ['./Schwartze_Paper1/figure_graphics/' monkey '_figure5_Dir_DOM_p'], '-dpng') 
% if monkey == 'P'
%     fig2svg('./figure_graphics/figure5b.svg', gcf)
% elseif monkey == 'Q'
%     fig2svg('./figure_graphics/figure5c.svg', gcf)
% end


vel_mod_ratio = (vel_Modulation_cor(All_signif_units) ./ vel_Modulation_init(All_signif_units));

vel_mod_init_ratio = (vel_Modulation_init1(All_signif_units,:) ./ vel_Modulation_init2(All_signif_units,:));
vel_mod_cor_ratio = (vel_Modulation_cor1(All_signif_units,:) ./ vel_Modulation_cor2(All_signif_units,:));

med_vel_mod_ratio = median(vel_mod_ratio(:));
disp(['Median Cor-Init FR ratio: ' num2str(med_vel_mod_ratio)])

init_prctile = prctile(vel_mod_init_ratio(:),[25,75]);
disp(['Init 25-75% ratio: ' num2str(init_prctile)])
cor_prctile = prctile(vel_mod_cor_ratio(:),[25,75]);
disp(['Cor 25-75% ratio: ' num2str(cor_prctile)])



%Mann-Whitney U test.
[p,h,stats] = ranksum(vel_mod_init_ratio(:),vel_mod_ratio(:));
[p,h,stats] = ranksum(vel_mod_cor_ratio(:),vel_mod_ratio(:));
%%
%Bar Graph below
cutoffs = 2.^(-2:0.5:2); % [-Inf,0.25,0.5,.75,1,1.5,2,4,Inf]; %[-Inf,0.5,1,2,Inf]
ratio_counts = histcounts(vel_mod_ratio,cutoffs);

% figure
% bar(ratio_counts, 'k');
% set(gca, 'XTickLabel', {'<0.5x', '0.5-1x', '1-2x', '>2x'})
% ylabel('Spiking Unit Count')
% set(gca, 'FontSize', 16)
% xlabel('Corrective to Initial FR ratio')
% title(['Monkey ' monkey])
% fig2svg('./figure_graphics/figure5f.svg', gcf)

ratio_init_counts = histcounts(vel_mod_init_ratio(:),cutoffs);
ratio_init_counts = sum(ratio_counts).*ratio_init_counts./sum(ratio_init_counts);
% figure
% bar(ratio_init_counts, 'FaceColor', init_color);

ratio_cor_counts = histcounts(vel_mod_cor_ratio(:),cutoffs);
ratio_cor_counts = sum(ratio_counts).*ratio_cor_counts./sum(ratio_cor_counts);
% figure
% bar(ratio_counts, 'k');
% hold on
% plot(1:length(ratio_init_counts),ratio_init_counts,'Color', init_color)
% plot(1:length(ratio_cor_counts),ratio_cor_counts,'Color', cor_color)
% bar(ratio_cor_counts, 'FaceColor', cor_color);

% figure
% plot(cumsum(ratio_counts))
% hold on
% plot(cumsum(ratio_init_counts))
% plot(cumsum(ratio_cor_counts))

% figure
% scatter(log2(sort(vel_mod_ratio)), 1:length(vel_mod_ratio))
% hold on
% plot(cutoffs(2:end), cumsum(ratio_init_counts))
% plot(cutoffs(2:end), cumsum(ratio_cor_counts))

n_signif_units = length(vel_mod_ratio);
n_bootstrap = length(vel_mod_init_ratio(:));


figure
plot(log2(sort(vel_mod_init_ratio(:))), (1:n_bootstrap)./1000,  'Color', init_color, 'LineWidth', 3)
hold on
plot(log2(sort(vel_mod_cor_ratio(:))), (1:n_bootstrap)./1000, 'Color', cor_color, 'LineWidth', 3)
scatter(log2(sort(vel_mod_ratio)), 1:length(vel_mod_ratio), 'MarkerEdgeColor', 'k')
xlim([-2,2])
ylim([0,n_signif_units])
line([1,1]*log2(med_vel_mod_ratio), [0,n_signif_units/2], 'Color', 'k', 'LineWidth', 2.5)
line([0,0], [0,n_signif_units/2], 'Color', 'k', 'LineWidth', 1.5, 'LineStyle','--')
scatter(log2(med_vel_mod_ratio), n_signif_units/2, 80, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w')
scatter(0, n_signif_units/2, 80, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
if monkey == 'P'
    text(log2(med_vel_mod_ratio)-0.1, -20, [num2str(med_vel_mod_ratio, '%1.2f'),'x'], 'FontSize', 14)
elseif monkey == 'Q'
    text(log2(med_vel_mod_ratio)-0.1, -16, [num2str(med_vel_mod_ratio, '%1.2f'),'x'], 'FontSize', 14)
end
    x_tick = -2:2;
x_tick_label = 2.^x_tick;
x_tick_label = {'0.25x', '0.5x', '1x', '2x', '4x'};
set(gca, 'XTick', x_tick)
set(gca, 'XTickLabel', x_tick_label)
xlabel('Corrective to Initial DOM ratio')
y_tick = linspace(0,n_signif_units,5);
y_tick_label = {'0', '25', '50', '75', '100'};
set(gca, 'YTick', y_tick)
set(gca, 'YTickLabel', y_tick_label)
ylabel('Percentile')
set(gca, 'FontSize', 14)
%fig2svg(['./figure_graphics/figure5f_' monkey '.svg'], gcf)
% if monkey == 'P'
%     fig2svg('./figure_graphics/figure5f.svg', gcf)
% elseif monkey == 'Q'
%     fig2svg('./figure_graphics/figure5g.svg', gcf)
% end
