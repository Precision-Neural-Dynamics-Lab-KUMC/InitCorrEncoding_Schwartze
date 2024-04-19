%% 
% Load data

clear variables
monkey = 'Q';

load(['\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results_StatsFinal'])


All_signif_units = squeeze(cat(3,signif_units{:}));
All_pref_dir_init = cat(2,pref_dir_init{:});
All_pref_dir_cor  = cat(2,pref_dir_cor{:});

for i = 1:12
    if i == 1
         All_pref_dir_init1= squeeze(pref_dir_init1{i});
    else
        All_pref_dir_init1 = [All_pref_dir_init1;squeeze(pref_dir_init1{i})];
    end
end
for i = 1:12
    if i == 1
         All_pref_dir_init2= squeeze(pref_dir_init2{i});
    else
        All_pref_dir_init2 = [All_pref_dir_init2;squeeze(pref_dir_init2{i})];
    end
end


for i = 1:12
    if i == 1
         All_pref_dir_cor1= squeeze(pref_dir_cor1{i});
    else
        All_pref_dir_cor1 = [All_pref_dir_cor1;squeeze(pref_dir_cor1{i})];
    end
end

for i = 1:12
    if i == 1
         All_pref_dir_cor2= squeeze(pref_dir_cor2{i});
    else
        All_pref_dir_cor2 = [All_pref_dir_cor2;squeeze(pref_dir_cor2{i})];
    end
end


PrefDirDiff=abs(wrapTo180(All_pref_dir_init-All_pref_dir_cor));


perc_neurons_lt_45deg = 100*(sum(PrefDirDiff(All_signif_units)<45)/sum(All_signif_units));

disp([num2str(100-perc_neurons_lt_45deg) '% Pref Dir > 45 degrees'])

%Bootstrapping
num_iter = 10000;
PrefDirData = PrefDirDiff(All_signif_units);
n_units = length(PrefDirData);
bootstrap_perc_neurons_lt_45deg = NaN(1,num_iter);
for k = 1:num_iter 
    rand_samp = randi(n_units, size(PrefDirData));
    bootstrap_perc_neurons_lt_45deg(k) = 100*(sum(PrefDirData(rand_samp)<45)/sum(All_signif_units));
    bootstrap_median_pd(k) = median(PrefDirData(rand_samp));

end
CI_perc_neurons_lt_45deg = prctile(bootstrap_perc_neurons_lt_45deg, [2.5,97.5]);
disp([num2str(100-CI_perc_neurons_lt_45deg(2)) '-' num2str(100-CI_perc_neurons_lt_45deg(1))  '% Pref Dir > 45 degrees'])

%% 42.4899% Pref Dir > 45 degrees
%% 38.9716-46.0081% Pref Dir > 45 degrees

CI_bootstrap_median_pd = prctile(bootstrap_median_pd, [2.5,97.5]);
disp([num2str(CI_bootstrap_median_pd(2)) '-' num2str(CI_bootstrap_median_pd(1))  '% Pref Dir > 45 degrees'])

disp(median(PrefDirDiff(All_signif_units)));

median_PrefDirDiff = median(PrefDirDiff(All_signif_units));


%Null Hypothesis 1: Neurons do not change PD at all
%initial submovements 
initPrefDirDiff = abs(wrapTo180(All_pref_dir_init1(All_signif_units,:)-All_pref_dir_init2(All_signif_units,:)));
median_initPrefDirDiff = median(initPrefDirDiff(:));

CI_median_initPrefDirDiff = prctile(initPrefDirDiff(:), [2.5,97.5]);
disp([num2str(CI_median_initPrefDirDiff(1)) '-' num2str(CI_median_initPrefDirDiff(2))  ' Degrees by 95% CI'])

median_median_initPrefDirDiff = median(median_initPrefDirDiff);

%corrective submovements
corPrefDirDiff = abs(wrapTo180(All_pref_dir_cor1(All_signif_units,:)-All_pref_dir_cor2(All_signif_units,:)));
median_corPrefDirDiff = median(corPrefDirDiff(:));

CI_median_corPrefDirDiff = prctile(corPrefDirDiff(:), [2.5,97.5]);
disp([num2str(CI_median_corPrefDirDiff(1)) '-' num2str(CI_median_corPrefDirDiff(2))  ' Degrees by 95% CI'])

median_median_corPrefDirDiff = median(median_corPrefDirDiff);



%Null Hypothesis 2: Randomizing corrective to match with initial submovements, all about the same
% Bootstrapping 
num_iter = 1000;
for k = 1:num_iter 
    
    rand_All_pref_dir_cor=randperm(length(All_pref_dir_cor));
    rand_PrefDirDiff = abs(wrapTo180(All_pref_dir_init-All_pref_dir_cor(rand_All_pref_dir_cor)));
    median_bootstrap_null(k,:) = median(rand_PrefDirDiff);

end

CI_median_bootstrap_null = prctile(median_bootstrap_null, [2.5,97.5]);
disp([num2str(CI_median_bootstrap_null(1)) '-' num2str(CI_median_bootstrap_null(2))  ' Degrees by 95% CI'])

median_median_bootstrap_null = median(median_bootstrap_null);


figure 
histogram(PrefDirDiff(All_signif_units), 0:15:180, 'FaceColor', [0.3,0.3,0.3])
hold on
ha = gca;
patch([CI_bootstrap_median_pd(1) CI_bootstrap_median_pd(1) CI_bootstrap_median_pd(2) CI_bootstrap_median_pd(2)], [min(ylim) max(ylim+20) max(ylim+20) min(ylim)], [0.8 0.8 0.8], 'edgecolor', 'none');
histogram(PrefDirDiff(All_signif_units), 0:15:180, 'FaceColor', [0.3,0.3,0.3])
xline(median_PrefDirDiff, 'LineStyle', '-', 'LineWidth', 2, 'Color', 'k')
xline([median_median_bootstrap_null], 'LineStyle', '--', 'LineWidth', 2, 'Color', '#7E2F8E');
xline([median_median_initPrefDirDiff], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'b');
xline([median_median_corPrefDirDiff], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r');
set(gca, 'XTick', 0:30:180)
xlabel('Preferred Direction Difference (degrees)')
ylabel('Number of units')
set(gca, 'FontSize', 14)
% set(gcf,'Position',[0 0 900 600]);
title(['Monkey ' monkey])
figName = 'Figure4C';
% print(figName, '-dtiff')
% print(figName, '-dpdf', '-painters' )
% fig2svg('./figure_graphics/figure4c.svg', gcf)
