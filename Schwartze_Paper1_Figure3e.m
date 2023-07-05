%% 
% Load data

clear variables
monkey = 'P';

date_strings = {'20170630', '20170712', '20170703', '20170713', '20170720', '20170731', '20170705', '20170706', '20170714', '20170717', '20170801', '20170802'}; 


addpath('./Violinplot-Matlab-master')


load(['\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results'])


All_signif_units = squeeze(cat(3,signif_units{:}));

All_best_RMSE_init = cat(2,best_RMSE_init{:});
All_best_RMSE_cor = cat(2,best_RMSE_cor{:});


%%
% Graphs

init_color = [0,0,0.7];
cor_color = [0.7,0,0];

[correlation_matrix,~,correlation_matrix_lower, correlation_matrix_upper] = corrcoef(All_best_RMSE_init(All_signif_units), All_best_RMSE_cor(All_signif_units));
correlation = correlation_matrix(1, 2);
disp(['Corr: ' num2str(correlation) ', (' num2str(correlation_matrix_lower(1,2)) '-' num2str(correlation_matrix_upper(1,2)) ') 95%CI'])


figure
scatter(All_best_RMSE_init(All_signif_units), All_best_RMSE_cor(All_signif_units), 'MarkerEdgeColor', 'k')
axis square
line([0,25],[0,25], 'LineStyle', '--', 'Color', 'k')
xlim([0,25])
ylim([0,25])
xlabel('Initial Model RMSE', 'Color', init_color)
ylabel('Corrective Model RMSE', 'Color', cor_color)
set(gca, 'FontSize', 14)
set(gcf,'Position',[0 0 600 400]);
hold on;
% fit = polyfit(All_best_RMSE_init(All_signif_units), All_best_RMSE_cor(All_signif_units), 1);
% y_fit = polyval(fit, All_best_RMSE_init(All_signif_units));
% plot(All_best_RMSE_init(All_signif_units), y_fit, 'color','cyan');
% print(gcf, ['./Schwartze_Paper1/figure_graphics/' monkey 'figure5_RMSE_init-cor_p'], '-dpng')

%paired t-test
[h,p,ci,stats] = ttest(All_best_RMSE_init(All_signif_units), All_best_RMSE_cor(All_signif_units));

% 
