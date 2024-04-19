%% 
% Load data

clear variables
monkey = 'P';
date_strings = {'20170630', '20170712', '20170703', '20170713', '20170720', '20170731', '20170705', '20170706', '20170714', '20170717', '20170801', '20170802'}; 


%%%

    data_path = file_locations('Data', '\Project_Data\20160504_COT_precision\data_extracted\monk_p\COT_SpikesCombined\');


d=1;
    
    load([data_path monkey '_Spikes_' date_strings{d} '-data.mat'])  %, 'TrialInfo', 'TrialSettings', 'speedPeaksTroughs', 'speedPeaksTroughs_i', 'speedPeaksRelProminence', 'speedPeaksPeakWidth', ...
%     'JoystickPos_disp','CursorSpeed', 'CursorSpeedRaw', 'TrialSpikeBins', 'TrialFiringRate', 'jPCsPhase', ...
%     'dataMask_peakVel', 'dataMask_zeroCross', 'trial_ids_peakVel', 'aligned_peaks', 'zeroCross_peakVel_i', 'initTrial_flag', 'smallInit_flag', 'withinTarg_peakVel_flag', ...
%     'b_regress_zeroCross', 'bint_regress_zeroCross', 'stats_regress_zeroCross', 'p_ttest_zeroCross', 'ci_ttest_zeroCross', 'stats_ttest_zeroCross')
%     


    All_TrialInfo(d) = TrialInfo;
    All_TrialSettings(d) = TrialSettings;

unique_trial_target = unique(All_TrialInfo(1).trial_target);



%%
center_target_x = All_TrialSettings(1).center_radii(1)*cos(0:(pi/100):(2*pi));
center_target_y = All_TrialSettings(1).center_radii(1)*sin(0:(pi/100):(2*pi));

figure
for k = 1:3
    subplot(1,3,k)
%     set(gcf,'CurrentAxes', ha{1}(k,1))
hold on
for targ = (k-1)*8+(1:8)
    wysiwyg

    draw_target(double(All_TrialSettings(1).target_angle(unique_trial_target(targ))), 2*pi/double(All_TrialSettings(1).num_display_targets(unique_trial_target(targ))), double(All_TrialSettings(1).inner_target_radii(unique_trial_target(targ))), double(All_TrialSettings(1).outer_target_radii(unique_trial_target(targ)))) 
    plot(center_target_x, center_target_y, 'Color', 'k')
    set(gcf,'color','w');
    axis square
    h=gcf;
    set(h,'Units','Points','Position',[0 0 600 400]); 

axis equal
axis off
end
end
% fig2svg(['./figure_graphics/figure1b.svg'], gcf)
