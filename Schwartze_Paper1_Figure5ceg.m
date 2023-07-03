%% 
% Load data

clear variables

monkey = 'Q';
date_strings = {'20180425', '20180426', '20180509', '20180510', '20180529', '20180530', '20180418', '20180419', '20180503', '20180507', '20180619', '20180620'};

load(['\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results'])


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


%Velocity regression coefficients depth of modulation
vel_Modulation_init = sqrt(sum(All_Coefficients_init(3:4,:).^2,1));
vel_Modulation_cor  = sqrt(sum(All_Coefficients_cor(3:4,:).^2,1));
    
%Data regression coefficients depth of modulation
data_Modulation_init = All_avg_speed_init.*vel_Modulation_init;
data_Modulation_cor  = All_avg_speed_cor.*vel_Modulation_cor;

%%
% Graphs

init_color = [0,0,0.7];
cor_color = [0.7,0,0];


% Calculate vel_mod correlation coefficient
[vel_mod_correlation_matrix,vmP,vmRL,vmRU] = corrcoef(vel_Modulation_init(All_signif_units), vel_Modulation_cor(All_signif_units));
vel_mod_correlation = vel_mod_correlation_matrix(1, 2);
disp(['Correlation init vs. corr speed modifed: ' num2str(vel_mod_correlation), ' 95% CI: ' num2str(vmRL(1,2)) '-' num2str(vmRU(1,2))])
%% Correlation init vs. corr speed modifed: 0.65515 95% CI: 0.60768-0.69796

figure
scatter(vel_Modulation_init(All_signif_units), vel_Modulation_cor(All_signif_units), 'MarkerEdgeColor', 'k')
hold on
line([0 7.5], [0 15], 'LineStyle', '--', 'Color', 'k')
line([0 15], [0 15], 'LineStyle', '--', 'Color', 'k')
line([0 15], [0 7.5], 'LineStyle', '--', 'Color', 'k')
xlim([-0.1,15])
ylim([-0.1,15])
text(13.2,7.9,'0.5x', 'FontSize', 12)
text(14.1,13.8,'1x', 'FontSize', 12)
text(7.5,14.6,'2x', 'FontSize', 12)
axis square
xlabel('Initial Velocity Depth of Mod.', 'Color', init_color)
ylabel('Corrective Velocity Depth of Mod.', 'Color', cor_color)
set(gca, 'FontSize', 14)
set(gcf,'Position',[0 0 600 400]);
title(['Monkey ' monkey])
% hold on;
% fit = polyfit(vel_Modulation_init(All_signif_units), vel_Modulation_cor(All_signif_units), 1);
% y_fit = polyval(fit, vel_Modulation_init(All_signif_units));
% plot(vel_Modulation_init(All_signif_units), y_fit, 'Color' , 'cyan');
% print(gcf, ['./Schwartze_Paper1/figure_graphics/' monkey '_figure5_Vel_DOM_p'], '-dpng')  




% % Calculate data_mod correlation coefficient
% data_mod_correlation_matrix = corrcoef(data_Modulation_init(All_signif_units), data_Modulation_cor(All_signif_units));
% data_mod_correlation = data_mod_correlation_matrix(1, 2);


figure
scatter(data_Modulation_init(All_signif_units), data_Modulation_cor(All_signif_units), 'MarkerEdgeColor', 'k')
hold on
line([0 15], [0 15], 'LineStyle', '--', 'Color', 'k')
xlim([-0.1,15])
ylim([-0.1,15])
axis square
xlabel('Initial Direction Depth of Mod.', 'Color', init_color)
ylabel('Corrective Direction Depth of Mod.', 'Color', cor_color)
set(gca, 'FontSize', 14)
set(gcf,'Position',[0 0 600 400]);
title(['Monkey ' monkey])
% hold on;
% fit = polyfit(data_Modulation_init(All_signif_units), data_Modulation_cor(All_signif_units), 1);
% y_fit = polyval(fit, data_Modulation_init(All_signif_units));
% plot(data_Modulation_init(All_signif_units), y_fit, 'Color','cyan');
% print(gcf, ['./Schwartze_Paper1/figure_graphics/' monkey '_figure5_Dir_DOM_p'], '-dpng') 
% 

vel_mod_ratio = (vel_Modulation_cor(All_signif_units) ./ vel_Modulation_init(All_signif_units));

%%
%Bar Graph below
ratio_counts = histcounts(vel_mod_ratio,[-Inf,0.5,1,2,Inf]);

figure
bar(ratio_counts, 'k');
set(gca, 'XTickLabel', {'<0.5x', '0.5-1x', '1-2x', '>2x'})
ylabel('Spiking Unit Count')
set(gca, 'FontSize', 16)
xlabel('Corrective to Initial FR ratio')
title(['Monkey ' monkey])
