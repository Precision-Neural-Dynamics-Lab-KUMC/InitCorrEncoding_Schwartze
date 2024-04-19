clear variables
% monkey = 'P';
% % date_strings = {'20170630'};
% date_strings = {'20170630', '20170712', '20170703', '20170713', '20170720', '20170731', '20170705', '20170706', '20170714', '20170717', '20170801', '20170802'};

monkey = 'Q';
date_strings = {'20180425', '20180426', '20180509', '20180510', '20180529', '20180530', '20180418', '20180419', '20180503', '20180507', '20180619', '20180620'};
% % Regression and data organization
run_analysis = true;
% run_analysis = false;



if run_analysis

    for d = 1:length(date_strings)
        if monkey == 'P'
            data_location  = 'C:\Users\kevin\Documents\MatLAb\Rouse\COT\monk_p\';
            if ~exist (data_location)
                data_location = 'C:\Users\legob\Documents\MATLAB\rouse\monk_p\';
            end
            if ~exist (data_location)
                data_location  = 'R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_extracted\monk_p\COT_SpikesCombined\';
            end
            if ~exist (data_location)
                data_location  = '\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_extracted\monk_p\COT_SpikesCombined\';
            end

        elseif monkey == 'Q'
            data_location  = 'C:\Users\kevin\Documents\MatLAb\Rouse\COT\monk_q\';
            if ~exist (data_location)
                data_location = 'C:\Users\legob\Documents\MATLAB\rouse\monk_q\';
            end
            if ~exist (data_location)
                data_location  = 'R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_extracted\monk_q\COT_SpikesCombined\';
            end
            if ~exist (data_location)
                data_location  = '\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_extracted\monk_q\COT_SpikesCombined\';
            end
        end

        load([data_location monkey '_Spikes_' date_strings{d}  '-data.mat'], 'CursorPos', 'TrialInfo', 'TrialSettings', 'CursorSpeedFilt', 'PeakInfo', 'PeakSettings');

        AllCursorPos{d} = CursorPos;
        %     PeakInit{d} = PeakInfo.speedPeaksTroughs_i(PeakInfo.initPeak_flag,2);
        %     day_index{d}  =  d*ones(size(PeakInit{d}));
        day_index{d}  =  d*ones(size(CursorPos,1));
        %     PeakInitTrials{d} = PeakInfo.trial_ids_peakVel(PeakInfo.initPeak_flag);
        %     PeakFirstCorr{d} = nan(size(PeakInit{d}));
        %     for tr = 1:length(PeakInitTrials{d})
        %         curr_peaks = PeakInfo.speedPeaksTroughs_i(PeakInfo.trial_ids_peakVel==PeakInitTrials{d}(tr),2);
        %         tmp = min(curr_peaks(curr_peaks>PeakInit{d}(tr)));
        %         if ~isempty(tmp)
        %             PeakFirstCorr{d}(tr) = tmp;
        %         end
        %     end


        %     data_before = 30;
        %     data_after = 60;
        %     initPeak_mask{d} = false(length(PeakInitTrials{d}),size(CursorSpeedRaw,2));
        %     for tr = 1:length(PeakInitTrials{d})
        %      initPeak_mask{d}(tr,PeakInit{d}(tr) + (-data_before:data_after)) = true;
        %     end
        AllCursorSpeed{d}=CursorSpeedFilt;
        AllTrialTargets{d} = TrialInfo.trial_target;
        %     CursorSpeed_initPeak{d} = getMaskedData(AllCursorSpeed{d}, initPeak_mask{d}, PeakInitTrials{d});
        %     CursorPos_initPeak{d} = getMaskedData(AllCursorPos{d}, initPeak_mask{d}, PeakInitTrials{d});

        %     AllPeakInfo{d} = PeakInfo;
        %     InitialTargets_initPeak{d} = AllTrialTargets{d}(PeakInitTrials{d},:);
        AllAlignSamples{d} = TrialInfo.align_samples;
        All_trial_ids_peakVel{d} = PeakInfo.trial_ids_peakVel;
        All_speedPeaksTroughs{d} = PeakInfo.speedPeaksTroughs;
        All_speedPeaksTroughs_i{d} = PeakInfo.speedPeaksTroughs_i;
        All_initPeak_flag{d} = PeakInfo.initPeak_flag;
        day_peak_index{d}  =  d*ones(size(PeakInfo.initPeak_flag,1));
    end

    load([data_location monkey '_Spikes_' date_strings{1}  '-data.mat'], 'TrialSettings')

    AllAlignSamples = cat(1,AllAlignSamples{:});
    num_trials_by_day = cellfun(@(x) size(x,1), AllCursorPos);
    total_trials_by_day = cumsum(num_trials_by_day);
    for d = 1:length(All_trial_ids_peakVel)
        if d==1
            All_trial_ids_peakVel{d} = All_trial_ids_peakVel{d};
        else
            All_trial_ids_peakVel{d} = total_trials_by_day(d-1) + All_trial_ids_peakVel{d};
        end
    end

    max_sample = max(cellfun(@(x) size(x,2), AllCursorPos));
    CursorPos_Trial = nan(total_trials_by_day(end), max_sample, 2);
    CursorSpeed_Trial = nan(total_trials_by_day(end), max_sample);
    for d = 1:length( AllCursorPos)
        if d==1
            CursorPos_Trial(1:total_trials_by_day(d), 1:size(AllCursorPos{d},2),:) = AllCursorPos{d};
            CursorSpeed_Trial(1:total_trials_by_day(d), 1:size(AllCursorSpeed{d},2)) = AllCursorSpeed{d};
        else
            CursorPos_Trial((total_trials_by_day(d-1)+1):total_trials_by_day(d), 1:size(AllCursorPos{d},2),:) = AllCursorPos{d};
            CursorSpeed_Trial((total_trials_by_day(d-1)+1):total_trials_by_day(d), 1:size(AllCursorPos{d},2)) = AllCursorSpeed{d};
        end
    end

    All_trial_ids_peakVel = cat(1,All_trial_ids_peakVel{:});
    All_initPeak_flag = cat(1,All_initPeak_flag{:});
    All_speedPeaksTroughs_i = cat(1, All_speedPeaksTroughs_i{:});
    AllTrialTargets = cat(1,AllTrialTargets{:});

    % time_vals = 10*(-data_before:data_after);
    % PeakInit = cat(1,PeakInit{:});
    % PeakInitTrials = cat(1,PeakInitTrials{:});
    % PeakFirstCorr  = cat(1,PeakFirstCorr{:});
    % CursorSpeed_initPeak  = cat(1,CursorSpeed_initPeak{:});
    % CursorPos_initPeak  = cat(1,CursorPos_initPeak{:});
    % InitialTargets_initPeak = cat(1,InitialTargets_initPeak{:});

end
EndIndex = find(strcmpi('End', TrialSettings.align_names));
sample_before = 30;
peak_sample_before = 10;
peak_sample_after = 10;
initPeakWholeTrial_mask = false(length(All_trial_ids_peakVel),size(CursorSpeed_Trial,2));
Peak_mask = false(length(All_trial_ids_peakVel),size(CursorSpeed_Trial,2));
for tr = 1:size(All_speedPeaksTroughs_i,1)
    if ~isnan(AllAlignSamples(All_trial_ids_peakVel(tr),EndIndex))
        initPeakWholeTrial_mask(tr,(All_speedPeaksTroughs_i(tr,2)-sample_before):(AllAlignSamples(All_trial_ids_peakVel(tr),EndIndex)-10)) = true;
    else
        %A few trials don't end, get all the data (these trials won't
        %actually be ploted
        initPeakWholeTrial_mask(tr,(All_speedPeaksTroughs_i(tr,2)-sample_before):end) = true;
    end
    Peak_mask(tr,All_speedPeaksTroughs_i(tr,2)+(-peak_sample_before:peak_sample_after)) = true;
end
CursorSpeed_Trial_peakVel = getMaskedData(CursorSpeed_Trial, initPeakWholeTrial_mask, All_trial_ids_peakVel);
CursorPos_Trial_peakVel = getMaskedData(CursorPos_Trial, initPeakWholeTrial_mask, All_trial_ids_peakVel);
CursorSpeed_peakVel = getMaskedData(CursorSpeed_Trial, Peak_mask, All_trial_ids_peakVel);
CursorPos_peakVel = getMaskedData(CursorPos_Trial, Peak_mask, All_trial_ids_peakVel);

center_target_x = TrialSettings(1).center_radii(1)*cos(0:(pi/100):(2*pi));
center_target_y = TrialSettings(1).center_radii(1)*sin(0:(pi/100):(2*pi));
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

targ = 1;
curr_trials = find(AllTrialTargets==(targ+16));
curr_trials = curr_trials(randi(length(curr_trials),5,1));

unique_trial_target = unique(AllTrialTargets);

figure('PaperPosition', [1,1,7.5, 3])
wysiwyg
ha = my_subplot([1,2,3],[1,1.25 1],[1,1,1],[0.25, .01, 0.25], [0.01,0.2,0.1], true);
ha{2}(1,1).Position(2) = 0.2;
ha{2}(1,1).Position(1) = 0.37;
rng(45)

for targ = 1:8 %:8
    set(gcf,'CurrentAxes', ha{1}(1,1))
    hold on
    curr_trials = find(AllTrialTargets==(targ+16));
    tmp_flag = false;
    while ~tmp_flag
        curr_trials = curr_trials(randi(length(curr_trials),6,1));
        if  max(max(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,50:end))) < max(max(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,1:50)))
            if all(all(isnan(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,400:end))))
                tmp_flag =true;
            end
        end
    end

    plot(CursorPos_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,:,1)', CursorPos_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,:,2)', 'Color', [0.7 0.7 0.7])

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag);
    plot(CursorPos_peakVel(curr_peaks, :,1)', CursorPos_peakVel(curr_peaks, :,2)', 'Color', init_color)
    scatter(CursorPos_peakVel(curr_peaks, peak_sample_before+1, 1), CursorPos_peakVel(curr_peaks, peak_sample_before+1, 2), 20, 's', 'MarkerFaceColor', 1-0.6*(1-init_color), 'MarkerEdgeColor', 1-0.6*(1-init_color))

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&~All_initPeak_flag);
    plot(CursorPos_peakVel(curr_peaks, :,1)', CursorPos_peakVel(curr_peaks, :,2)', 'Color', cor_color)
    scatter(CursorPos_peakVel(curr_peaks, peak_sample_before+1, 1), CursorPos_peakVel(curr_peaks, peak_sample_before+1, 2), 20, 's', 'MarkerFaceColor', 1-0.6*(1-cor_color), 'MarkerEdgeColor', 1-0.6*(1-cor_color))


    draw_target(double(TrialSettings.target_angle(unique_trial_target(targ+16))), 2*pi/8, double(TrialSettings.inner_target_radii(unique_trial_target(targ+16))), double(TrialSettings.outer_target_radii(unique_trial_target(targ+16))))

    if targ ==8
        plot(center_target_x, center_target_y, 'Color', 'k')
        axis square
        axis equal
        axis off
    end


    set(gcf,'CurrentAxes', ha{2}(1,1))
    hold on
    plot(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,:,1)', 'Color',  [0.2 0.2 0.2])

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag);
    for p = 1:length(curr_peaks)
        curr_t =All_speedPeaksTroughs_i(curr_peaks(p),2) - All_speedPeaksTroughs_i( All_trial_ids_peakVel==All_trial_ids_peakVel(curr_peaks(p))&All_initPeak_flag,2 ) + sample_before +1 ;
        scatter(curr_t, CursorSpeed_peakVel(curr_peaks(p), peak_sample_before+1), 's', 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
    end

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&~All_initPeak_flag);
    for p = 1:length(curr_peaks)
        curr_t =All_speedPeaksTroughs_i(curr_peaks(p),2) - All_speedPeaksTroughs_i( All_trial_ids_peakVel==All_trial_ids_peakVel(curr_peaks(p))&All_initPeak_flag,2 ) + sample_before + 1;
        scatter(curr_t, CursorSpeed_peakVel(curr_peaks(p), peak_sample_before+1), 's', 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
    end
end
xlim([0, 320])
line(get(gca,'XLim'), [1 1] *250, 'LineStyle', '--', 'Color', 'k')
ylabel( 'Cursor speed, pixels/s')
set(gca, 'XTick', 30+(0:100:200))
set(gca, 'FontSize', 14)
xlabel('Time, ms')

x_ticks = get(gca, 'XTickLabel');
x_ticks{1} = '0';
x_ticks{2} = '1000';
x_ticks{3} = '2000';
% x_ticks{4} = '3000';
set(gca, 'XTickLabel', x_ticks)

annotation('textbox',[0.01 0.94 0.05 0.03], 'String', 'A', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')
annotation('textbox',[0.26 0.94 0.05 0.03], 'String', 'B', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')
annotation('textbox',[0.70 0.94 0.05 0.03], 'String', 'C', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')

set(gcf,'CurrentAxes', ha{3}(1,1))
axis square
axis equal
axis off

figName = [monkey '_Figure2ab'];
% print(figName, '-dtiff')
% print(figName, '-dpdf', '-vector' )

%%Figure 2C is plotted separately so you can save as raster graphic as too
%%many data points for vector

figure('PaperPosition', [1,1,7.5, 3])
wysiwyg
ha = my_subplot([1,2,3],[1,1.25 1],[1,1,1],[0.25, .01, 0.25], [0.01,0.2,0.1], true);
ha{2}(1,1).Position(2) = 0.2;
ha{2}(1,1).Position(1) = 0.37;
rng(45)
set(gcf,'CurrentAxes', ha{1}(1,1))
axis square
axis equal
axis off
set(gcf,'CurrentAxes', ha{2}(1,1))
axis off
set(gcf,'CurrentAxes', ha{3}(1,1))
submovementsToPlot = false(size(All_initPeak_flag));
submovementsToPlot(1:10:12637) = true;
submovementsToPlotCor = ~All_initPeak_flag & submovementsToPlot;
plot(CursorPos_peakVel(submovementsToPlotCor,[1,end],1)',CursorPos_peakVel(submovementsToPlotCor,[1,end],2)','Color',cor_color)
hold on
scatter(CursorPos_peakVel(submovementsToPlotCor,1,1)',CursorPos_peakVel(submovementsToPlotCor,1,2)',5,'MarkerFaceColor',cor_color,'MarkerEdgeColor',[0,0,0])

submovementsToPlotInit = All_initPeak_flag & submovementsToPlot;
plot(CursorPos_peakVel(submovementsToPlotInit,[1,end],1)',CursorPos_peakVel(submovementsToPlotInit,[1,end],2)','Color',init_color)
scatter(CursorPos_peakVel(submovementsToPlotInit,1,1)',CursorPos_peakVel(submovementsToPlotInit,1,2)',5,'MarkerFaceColor',init_color,'MarkerEdgeColor',[0,0,0])
axis square
axis equal
axis off

annotation('textbox',[0.01 0.94 0.05 0.03], 'String', 'A', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')
annotation('textbox',[0.26 0.94 0.05 0.03], 'String', 'B', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')
annotation('textbox',[0.70 0.94 0.05 0.03], 'String', 'C', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')

figName = [monkey '_Figure2c'];
print(figName, '-dtiff')
print(figName, '-dpdf', '-vector' )


figure('PaperPosition', [1,1,7.5, 3])
wysiwyg
ha = my_subplot([1,2,3],[1,1.25 1],[1,1,1],[0.25, .01, 0.25], [0.01,0.2,0.1], true);
ha{2}(1,1).Position(2) = 0.2;
ha{2}(1,1).Position(1) = 0.37;
rng(45)

for targ = 1:8 %:8
    set(gcf,'CurrentAxes', ha{1}(1,1))
    hold on
    curr_trials = find(AllTrialTargets==(targ+16));
    tmp_flag = false;
    while ~tmp_flag
        curr_trials = curr_trials(randi(length(curr_trials),6,1));
        if  max(max(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,50:end))) < max(max(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,1:50)))
            if all(all(isnan(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,400:end))))
                tmp_flag =true;
            end
        end
    end

    plot(CursorPos_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,:,1)', CursorPos_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,:,2)', 'Color', [0.7 0.7 0.7])

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag);
    plot(CursorPos_peakVel(curr_peaks, :,1)', CursorPos_peakVel(curr_peaks, :,2)', 'Color', init_color)
    scatter(CursorPos_peakVel(curr_peaks, peak_sample_before+1, 1), CursorPos_peakVel(curr_peaks, peak_sample_before+1, 2), 20, 's', 'MarkerFaceColor', 1-0.6*(1-init_color), 'MarkerEdgeColor', 1-0.6*(1-init_color))

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&~All_initPeak_flag);
    plot(CursorPos_peakVel(curr_peaks, :,1)', CursorPos_peakVel(curr_peaks, :,2)', 'Color', cor_color)
    scatter(CursorPos_peakVel(curr_peaks, peak_sample_before+1, 1), CursorPos_peakVel(curr_peaks, peak_sample_before+1, 2), 20, 's', 'MarkerFaceColor', 1-0.6*(1-cor_color), 'MarkerEdgeColor', 1-0.6*(1-cor_color))


    draw_target(double(TrialSettings.target_angle(unique_trial_target(targ+16))), 2*pi/8, double(TrialSettings.inner_target_radii(unique_trial_target(targ+16))), double(TrialSettings.outer_target_radii(unique_trial_target(targ+16))))

    if targ ==8
        plot(center_target_x, center_target_y, 'Color', 'k')
        axis square
        axis equal
        axis off
    end


    set(gcf,'CurrentAxes', ha{2}(1,1))
    hold on

    plot(CursorSpeed_Trial_peakVel(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag,:,1)', 'Color',  [0.2 0.2 0.2])

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&All_initPeak_flag);
    for p = 1:length(curr_peaks)
        curr_t =All_speedPeaksTroughs_i(curr_peaks(p),2) - All_speedPeaksTroughs_i( All_trial_ids_peakVel==All_trial_ids_peakVel(curr_peaks(p))&All_initPeak_flag,2 ) + sample_before +1 ;
        scatter(curr_t, CursorSpeed_peakVel(curr_peaks(p), peak_sample_before+1), 's', 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
    end

    curr_peaks = find(ismember(All_trial_ids_peakVel,curr_trials)&~All_initPeak_flag);
    for p = 1:length(curr_peaks)
        curr_t =All_speedPeaksTroughs_i(curr_peaks(p),2) - All_speedPeaksTroughs_i( All_trial_ids_peakVel==All_trial_ids_peakVel(curr_peaks(p))&All_initPeak_flag,2 ) + sample_before + 1;
        scatter(curr_t, CursorSpeed_peakVel(curr_peaks(p), peak_sample_before+1), 's', 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
    end
end
xlim([0, 320])
line(get(gca,'XLim'), [1 1] *250, 'LineStyle', '--', 'Color', 'k')
ylabel( 'Cursor speed, pixels/s')
set(gca, 'XTick', 30+(0:100:200))
set(gca, 'FontSize', 14)
xlabel('Time, ms')

x_ticks = get(gca, 'XTickLabel');
x_ticks{1} = '0';
x_ticks{2} = '1000';
x_ticks{3} = '2000';
% x_ticks{4} = '3000';
set(gca, 'XTickLabel', x_ticks)

set(gcf,'CurrentAxes', ha{3}(1,1))
submovementsToPlot = false(size(All_initPeak_flag));
submovementsToPlot(1:10:12637) = true;
submovementsToPlotCor = ~All_initPeak_flag & submovementsToPlot;
plot(CursorPos_peakVel(submovementsToPlotCor,[1,end],1)',CursorPos_peakVel(submovementsToPlotCor,[1,end],2)','Color',cor_color)
hold on
scatter(CursorPos_peakVel(submovementsToPlotCor,1,1)',CursorPos_peakVel(submovementsToPlotCor,1,2)',5,'MarkerFaceColor',cor_color,'MarkerEdgeColor',[0,0,0])

submovementsToPlotInit = All_initPeak_flag & submovementsToPlot;
plot(CursorPos_peakVel(submovementsToPlotInit,[1,end],1)',CursorPos_peakVel(submovementsToPlotInit,[1,end],2)','Color',init_color)
scatter(CursorPos_peakVel(submovementsToPlotInit,1,1)',CursorPos_peakVel(submovementsToPlotInit,1,2)',5,'MarkerFaceColor',init_color,'MarkerEdgeColor',[0,0,0])
axis square
axis equal
axis off
annotation('textbox',[0.01 0.94 0.05 0.03], 'String', 'A', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')
annotation('textbox',[0.26 0.94 0.05 0.03], 'String', 'B', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')
annotation('textbox',[0.70 0.94 0.05 0.03], 'String', 'C', 'EdgeColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center')
figName = [monkey '_Figure2'];
print(figName, '-dtiff')
print(figName, '-dpdf', '-vector' )



