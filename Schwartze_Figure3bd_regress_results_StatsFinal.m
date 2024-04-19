%% 
% Load data

clear variables

monkey = 'Q';
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
All_pref_dir_init = cat(2,pref_dir_init{:});
All_pref_dir_cor  = cat(2,pref_dir_cor{:});
All_best_R2_tot_init = squeeze(cat(3,best_R2_tot_init{:}));
All_best_R2_tot_cor = squeeze(cat(3,best_R2_tot_cor{:}));
All_best_R2_partial_init = cat(2,best_R2_partial_init{:});
All_best_R2_partial_cor  = cat(2,best_R2_partial_cor{:});



%%
% Graphs

init_color = [0,0,0.7];
cor_color = [0.7,0,0];

vel_color = [255,146,3]/255;
pos_color = [2,194,34]/255;
speed_color = [245,117,232]/255;
tot_color = [77,77,77]/255;


opp_vel_color = [77,58,250]/255;
opp_pos_color = [189,98,96]/255;
opp_speed_color = [232,245,117]/255;

%Selected using https://www.sessions.edu/color-calculator/
vel_color2 = [255,238,0]/255;
vel_color3 = [255,43,0]/255;
vel_colormap = create_color_map3(vel_color, vel_color2, vel_color3);
pos_color2 = [2,21,199]/255;
pos_color3 = [170,196,0]/255;
pos_colormap = create_color_map3(pos_color, pos_color2, pos_color3);
num_mag_levels = 20;
vel_colormap2d = zeros(size(vel_colormap,1),num_mag_levels,3);
for k = 1:size(vel_colormap2d,1)
    for c = 1:3
    vel_colormap2d(k,:,c) = linspace(vel_colormap(k,c),0.9,num_mag_levels);
    end
end

pos_colormap2d = zeros(size(pos_colormap,1),num_mag_levels,3);
for k = 1:size(pos_colormap2d,1)
    for c = 1:3
    pos_colormap2d(k,:,c) = linspace(pos_colormap(k,c),0.9,num_mag_levels);
    end
end
for c = 1:3
speed_colormap(:,c) = linspace(speed_color(c), 0.9, num_mag_levels);
end




plot_R2_partial_init = cat(1, permute(All_best_R2_tot_init(All_signif_units),[2,1]),  All_best_R2_partial_init(:,All_signif_units));
plot_R2_partial_init = plot_R2_partial_init(:);
plot_R2_partial_cor = cat(1, permute(All_best_R2_tot_cor(All_signif_units),[2,1]),  All_best_R2_partial_cor(:,All_signif_units));
plot_R2_partial_cor = plot_R2_partial_cor(:);
labels_R2_partial = repmat({'Total'; 'Position'; 'Velocity'; 'Speed'},[1,sum(All_signif_units)]);
labels_R2_partial = labels_R2_partial(:);
mean_R2_partial_init = mean(All_best_R2_partial_init(:,All_signif_units),2);
mean_R2_partial_cor  = mean(All_best_R2_partial_cor(:,All_signif_units),2);


%Plot partial R^2 for initial movements
figure('PaperPosition', 2*[1,1,2.5,2.5*0.75])  %Made at 2x size
wysiwyg;
vh = violinplot(plot_R2_partial_init, labels_R2_partial, 'GroupOrder', {'Total', 'Velocity', 'Position', 'Speed'});
vh(1).ViolinColor = {tot_color};
vh(2).ViolinColor = {vel_color};
vh(3).ViolinColor = {pos_color};
vh(4).ViolinColor = {speed_color};
ylabel('R^2')
set(gca, 'FontSize', 14)
ylim([0, 0.8])
% title(['Monkey ' monkey ' Initial Submovements'])
ax = axes('Position', [.6,.4,.4,.4]);
pievelinit=round((mean_R2_partial_init(2)/sum(mean_R2_partial_init)*100),1);
pieposinit=round((mean_R2_partial_init(1)/sum(mean_R2_partial_init)*100),1);
piespinit=round((mean_R2_partial_init(3)/sum(mean_R2_partial_init)*100),1);
ph = pie(mean_R2_partial_init([2,1,3]), {['Vel: ' num2str(pievelinit) '%'],[ 'Pos: ' num2str(pieposinit) '%'],[ 'Speed: ' num2str(piespinit) '%']});
ph(2).FontSize = 14;  ph(4).FontSize = 14;  ph(6).FontSize = 14;
%Pi labels manually adjusted
colormap([vel_color; pos_color; speed_color])
ph(4).Position(2) = ph(4).Position(2)+0.05;
ph(6).Position(2) = ph(6).Position(2)-0.05;
% set(gcf,'Position',[0 0 600 400]);
% print(gcf, ['./Schwartze_Paper1/figure_graphics/' monkey '_figure3_violin_init_p'], '-dpng')
% fig2svg(['./figure_graphics/' monkey '_figure3_violin_init_p_test2.svg'], gcf)

% %Plot partial R^2 for corrective movements
figure('PaperPosition', 2*[1,1,2.5,2.5*0.75])  %Made at 2x size
wysiwyg;
vh = violinplot(plot_R2_partial_cor, labels_R2_partial, 'GroupOrder', {'Total', 'Velocity', 'Position', 'Speed'});
vh(1).ViolinColor = {tot_color};
vh(2).ViolinColor = {vel_color};
vh(3).ViolinColor = {pos_color};
vh(4).ViolinColor = {speed_color};
ylabel('R^2')
% title(['Monkey ' monkey ' Corrective Submovements'])
set(gca, 'FontSize', 14)
ylim([0, 0.8])
ax = axes('Position', [.6,.4,.4,.4]);
pievelcor=round((mean_R2_partial_cor(2)/sum(mean_R2_partial_cor)*100),1);
pieposcor=round((mean_R2_partial_cor(1)/sum(mean_R2_partial_cor)*100),1);
piespcor=round((mean_R2_partial_cor(3)/sum(mean_R2_partial_cor)*100),1);
ph = pie(mean_R2_partial_cor([2,1,3]), {['Vel: ' num2str(pievelcor) '%'],[ 'Pos: ' num2str(pieposcor) '%'],[ 'Speed: ' num2str(piespcor) '%']});
ph(2).FontSize = 14;  ph(4).FontSize = 14;  ph(6).FontSize = 14;
%Pi labels manually adjusted
colormap([vel_color; pos_color; speed_color])
% set(gcf,'Position',[0 0 600 400]);
% print(gcf, ['./Schwartze_Paper1/figure_graphics/' monkey '_figure3_violin_cor_p'], '-dpng')
% fig2svg(['./figure_graphics/' monkey '_figure3_violin_cor_p_test2.svg'], gcf)


