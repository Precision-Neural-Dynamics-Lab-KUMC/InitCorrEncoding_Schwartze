%%
% Load data

clear variables

monkey = 'P';
date_strings = {'20170630', '20170712', '20170703', '20170713', '20170720', '20170731', '20170705', '20170706', '20170714', '20170717', '20170801', '20170802'};

% monkey = 'Q';
% date_strings = {'20180425', '20180426', '20180509', '20180510', '20180529', '20180530', '20180418', '20180419', '20180503', '20180507', '20180619', '20180620'};

% Regression and data organization
% run_analysis = true;
run_analysis = false;
% 
% addpath('./Violinplot-Matlab-master')

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

            M1_Arrays = {'G', 'H', 'I', 'J', 'K', 'L'};
            S1_Arrays = {'M','N', 'O', 'P'};
            PMv_Arrays = {'A+', 'A', 'C'};  %A+, A in PMv, EFXYZ in PMd (C and D PM in between)
            PMd_Arrays = {'E', 'F', 'X', 'Y', 'Z'};
            PRR_Arrays = {'U', 'V', 'V''', 'W'};
            PMo_Arrays = {'D'};
        elseif monkey == 'Q'
            data_location  = 'C:\Users\kevin\Documents\MatLAb\Rouse\COT\monk_q\';
            if ~exist (data_location)
                data_location = 'C:\Users\legob\Documents\MATLAB\rouse\monk_q\';
            end
            if ~exist (data_location)
                data_location  = 'R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_extracted\monk_q\COT_SpikesCombined\';
            end
            % load([data_location 'Q_Spikes_20180419-data.mat'])
            M1_Arrays = { 'H', 'I', 'J', 'K', 'alpha'};
            S1_Arrays = {'M','N', 'O', 'P'};
            PMv_Arrays = {'A', 'B', 'C', 'D'};  %ABCD in PMv, EFYZ in PMd
            PMd_Arrays = {'E', 'F', 'Y', 'Z'};
            PRR_Arrays = {'U', 'V', 'W', 'X'};
        end


        %%

        load([data_location monkey '_Spikes_' date_strings{d}  '-data.mat'])  %



        %Use these to pull out the area for a given spiking unit in
        %SpikeFiringRates or SpikeTimes
        M1_units_i = ismember(SpikeSettings.array_by_channel, M1_Arrays);
        S1_units_i = ismember(SpikeSettings.array_by_channel, S1_Arrays);
        PMv_units_i = ismember(SpikeSettings.array_by_channel, PMv_Arrays);
        PMd_units_i = ismember(SpikeSettings.array_by_channel, PMd_Arrays);
        PRR_units_i = ismember(SpikeSettings.array_by_channel, PRR_Arrays);

        %%Pull out firing rates in a time window around movement speed peaks
        start_samples_peakVel = -50;  %500 ms before movement peak speeds
        end_samples_peakVel = 30;  %300 ms after movement peak speeds
        mid_samples_peakVel = -start_samples_peakVel+1;
        time_array = 1000*(start_samples_peakVel:end_samples_peakVel)/SpikeSettings.samp_rate; %Time array in ms

        %PeakInfo.speedPeaksTroughs_i(:,2) is the sample of the maximum peak speed
        %   speedPeaksTroughs_i(:,1) is the trough before
        %   speedPeaksTroughs_i(:,3) is the trough after
        %PeakInfo.speedPeaksTroughs has the speed values
        dataMask_peakVel = zeros(size(PeakInfo.speedPeaksTroughs_i,1),size(SpikeFiringRates,2));
        for n = 1:size(PeakInfo.speedPeaksTroughs_i,1)
            curr_indexes = PeakInfo.speedPeaksTroughs_i(n,2) + (start_samples_peakVel:end_samples_peakVel);
            curr_indexes = curr_indexes(curr_indexes>0);
            dataMask_peakVel(n,curr_indexes) = 1;
        end
        Joystick_position = getMaskedData(JoystickPos_disp, dataMask_peakVel, PeakInfo.trial_ids_peakVel);

        start_pos = squeeze(Joystick_position(:,41,:)); %100 ms before peak vel
        end_pos = squeeze(Joystick_position (:,61,:)); %100 ms after peak vel
        timestep = 0.01; %sec

        [X,SD] = normalization(start_pos, end_pos, timestep);

        %All firing rates from 500ms before movement peak speeds until 300ms after
        SpikeFiringRates_peakVel = getMaskedData(SpikeFiringRates, dataMask_peakVel, PeakInfo.trial_ids_peakVel);

        % Select for brain region of study and time lag for regressions, if desired
        SpikeFiringRates_peakVel = SpikeFiringRates_peakVel(:,:,M1_units_i);

        %select firing rates for initial or corrective peaks only
        SpikeFiringRates_peakVel_init = SpikeFiringRates_peakVel(PeakInfo.initPeak_flag,:,:);
        SpikeFiringRates_peakVel_cor = SpikeFiringRates_peakVel(~PeakInfo.initPeak_flag,:,:);

        X_vel_angle = atan2d(X(:,4),X(:,3));
        X_pos_angle = atan2d(X(:,2),X(:,1));

        %select X values from normalization for inital and corrective peaks only
        X_init = X(PeakInfo.initPeak_flag,:);
        X_cor = X(~PeakInfo.initPeak_flag,:);
        X_vel_angle_init = X_vel_angle(PeakInfo.initPeak_flag,:);
        X_vel_angle_cor = X_vel_angle(~PeakInfo.initPeak_flag,:);
        X_pos_angle_init = X_pos_angle(PeakInfo.initPeak_flag,:);
        X_pos_angle_cor = X_pos_angle(~PeakInfo.initPeak_flag,:);


        %all initial and corrective peaks ran together, per neuron
        R2_partial = NaN(3,41,size(SpikeFiringRates_peakVel,3));
        R2_tot = NaN(1,41,size(SpikeFiringRates_peakVel,3));
        Coefficients = NaN(6,41,size(SpikeFiringRates_peakVel,3));
        p_vals = NaN(6,41,size(SpikeFiringRates_peakVel,3));
        for n = 1:size(SpikeFiringRates_peakVel,3)
            for t = 1:41
                [R2_partial(:,t,n), R2_tot(:,t,n), Coefficients(:,t,n), p_vals(:,t,n), curr_RMSE(t,n)] = Do_mutiregress(X,SpikeFiringRates_peakVel(:,(20+t),n));
            end
        end
        disp(['Day ' num2str(d) ' all regression done'])
        %
        %initial peaks only
        R2_partial_init = NaN(3,41,size(SpikeFiringRates_peakVel,3));
        R2_tot_init = NaN(1,41,size(SpikeFiringRates_peakVel,3));
        Coefficients_init = NaN(6,41,size(SpikeFiringRates_peakVel,3));
        p_vals_init = NaN(6,41,size(SpikeFiringRates_peakVel,3));
        for n = 1:size(SpikeFiringRates_peakVel,3)
            parfor t = 1:41
                [R2_partial_init(:,t,n), R2_tot_init(:,t,n), Coefficients_init(:,t,n), p_vals_init(:,t,n), curr_RMSE_init(t,n)] = Do_mutiregress(X_init,SpikeFiringRates_peakVel_init(:,(20+t),n));
            end
        end
        disp(['Day ' num2str(d) ' init regression done'])
        %
        %corretive peaks only
        R2_partial_cor = NaN(3,41,size(SpikeFiringRates_peakVel,3));
        R2_tot_cor = NaN(1,41,size(SpikeFiringRates_peakVel,3));
        Coefficients_cor = NaN(6,41,size(SpikeFiringRates_peakVel,3));
        p_vals_cor = NaN(6,41,size(SpikeFiringRates_peakVel,3));
        for n = 1:size(SpikeFiringRates_peakVel,3)
            parfor t = 1:41
                [R2_partial_cor(:,t,n), R2_tot_cor(:,t,n), Coefficients_cor(:,t,n), p_vals_cor(:,t,n), curr_RMSE_cor(t,n)] = Do_mutiregress(X_cor,SpikeFiringRates_peakVel_cor(:,(20+t),n));
            end
        end
        disp(['Day ' num2str(d) ' cor regression done'])

        All_X_init{d} = X_init;
        All_X_cor{d} = X_cor;
        All_SpikeFiringRates_peakVel_init{d} = SpikeFiringRates_peakVel_init;
        All_SpikeFiringRates_peakVel_cor{d} = SpikeFiringRates_peakVel_cor;
        All_Coefficients_init{d} = Coefficients_init;
        All_Coefficients_cor{d} = Coefficients_cor;

        RMSE{d} = curr_RMSE;
        RMSE_init{d} =  curr_RMSE_init;
        RMSE_cor{d} = curr_RMSE_cor;
        %
        std_FR_init = NaN(41,size(SpikeFiringRates_peakVel,3));
        std_FR_cor = NaN(41,size(SpikeFiringRates_peakVel,3));
        for n = 1:size(SpikeFiringRates_peakVel,3)
            for t = 1:41
                std_FR_init(t,n) = std(SpikeFiringRates_peakVel_init(:,(20+t),n));
                std_FR_cor(t,n)  = std(SpikeFiringRates_peakVel_cor(:,(20+t),n));
            end
        end

        avg_speed_init(d) = mean(sqrt(sum(X_init(:,3:4).^2,2)));
        avg_speed_cor(d) = mean(sqrt(sum(X_cor(:,3:4).^2,2)));

        [best_R2_tot_init{d}, best_init_t{d}] = max(R2_tot_init,[],2);
        [best_R2_tot_cor{d}, best_cor_t{d}] = max(R2_tot_cor,[],2);
        signif_units{d} = best_R2_tot_init{d}>.1 | best_R2_tot_cor{d}>.1;

        for n = 1:length(best_init_t{d})
            % best_SpikeFiringRates_peakVel_init(:,n) = SpikeFiringRates_peakVel_init(:,(20+best_init_t(n)),n);
            % best_SpikeFiringRates_peakVel_cor(:,n) = SpikeFiringRates_peakVel_cor(:,(20+best_cor_t(n)),n);
            best_Coefficients_init{d}(:,n) = Coefficients_init(:,best_init_t{d}(n),n);
            best_Coefficients_cor{d}(:,n) = Coefficients_cor(:,best_cor_t{d}(n),n);
            best_std_FR_init{d}(n) = std_FR_init(best_init_t{d}(n),n);
            best_std_FR_cor{d}(n)  = std_FR_init(best_cor_t{d}(n),n);
            pref_dir_init{d}(n) = atan2d(best_Coefficients_init{d}(5,n), best_Coefficients_init{d}(4,n));
            pref_dir_cor{d}(n) = atan2d(best_Coefficients_cor{d}(5,n), best_Coefficients_cor{d}(4,n));
            best_R2_partial_init{d}(:,n) = squeeze(R2_partial_init(:,best_init_t{d}(n),n));
            best_R2_partial_cor{d}(:,n)  = squeeze(R2_partial_cor(:,best_cor_t{d}(n),n));
            best_RMSE_init{d}(n) = RMSE_init{d}(best_init_t{d}(n),n);
            best_RMSE_cor{d}(n) = RMSE_cor{d}(best_cor_t{d}(n),n);
        end

        clearvars -except monkey date_strings d avg_speed_init avg_speed_cor best_R2_tot_init best_init_t best_R2_tot_cor best_cor_t signif_units best_Coefficients_init best_Coefficients_cor ...
            best_std_FR_init best_std_FR_cor pref_dir_init pref_dir_cor best_R2_partial_init best_R2_partial_cor RMSE RMSE_init RMSE_cor best_RMSE_init best_RMSE_cor All_X_init All_X_cor ...
            All_SpikeFiringRates_peakVel_init All_SpikeFiringRates_peakVel_cor All_Coefficients_init All_Coefficients_cor

        d
    end
%     save(['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results'])
    save(['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_all_data'])
else
%     load(['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results'])
    load(['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_all_data'])
    load(['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results_PDStats3'], 'best_Coefficients_init1')
end
% Q_regress_results_PDStats3
All_day_index = cellfun(@(x) ones(size(x,3),1), signif_units, 'UniformOutput', false);
for d = 1:length(signif_units)
    All_day_index{d} = d*ones(size(signif_units{d},3),1);
end
All_day_index = cat(1,All_day_index{:});


All_avg_speed_init = avg_speed_init(All_day_index);
All_avg_speed_cor  = avg_speed_cor(All_day_index);

All_signif_units = squeeze(cat(3,signif_units{:}));
All_best_Coefficients_init = cat(2,best_Coefficients_init{:});
Bootstrap_best_Coefficients_init = squeeze(cat(3,best_Coefficients_init1{:}));
All_best_Coefficients_cor = cat(2,best_Coefficients_cor{:});
All_pref_dir_init = cat(2,pref_dir_init{:});
All_pref_dir_cor  = cat(2,pref_dir_cor{:});
All_best_R2_tot_init = squeeze(cat(3,best_R2_tot_init{:}));
All_best_R2_tot_cor = squeeze(cat(3,best_R2_tot_cor{:}));
All_best_R2_partial_init = cat(2,best_R2_partial_init{:});
All_best_R2_partial_cor  = cat(2,best_R2_partial_cor{:});
All_best_RMSE_init = cat(2,best_RMSE_init{:});
All_best_RMSE_cor = cat(2,best_RMSE_cor{:});
All_Coefficients_init = cat(3,All_Coefficients_init{:});
All_Coefficients_cor = cat(3,All_Coefficients_cor{:});

%Velocity regression coefficients depth of modulation
vel_Modulation_init = sqrt(sum(All_best_Coefficients_init(3:4,:).^2,1));
vel_Modulation_cor  = sqrt(sum(All_best_Coefficients_cor(3:4,:).^2,1));

%Data regression coefficients depth of modulation
data_Modulation_init = All_avg_speed_init.*vel_Modulation_init;
data_Modulation_cor  = All_avg_speed_cor.*vel_Modulation_cor;


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





PrefDirDiff=abs(wrapTo180(All_pref_dir_init-All_pref_dir_cor));









mean_init_speed = mean(All_avg_speed_init);
mean_cor_speed = mean(All_avg_speed_cor);

fit_vectors_init = [ones(1,8); zeros(2,8); mean_init_speed*cosd(0:45:315); mean_init_speed*sind(0:45:315); mean_init_speed*ones(1,8)]';
fit_vectors_cor = [ones(1,8); zeros(2,8);  mean_cor_speed*cosd(0:45:315);  mean_cor_speed*sind(0:45:315); mean_cor_speed*ones(1,8)]';


fit_FR_init = pagemtimes(fit_vectors_init, All_Coefficients_init );
fit_FR_cor = pagemtimes(fit_vectors_cor, All_Coefficients_cor );


t_pt_for_pca = 16;  %150 ms before movement
fit_init_pca_coeff = pca(reshape(fit_FR_init(:,t_pt_for_pca ,All_signif_units), [], sum(All_signif_units)), 'Centered', true);
fit_cor_pca_coeff = pca(reshape(fit_FR_cor(:,t_pt_for_pca ,All_signif_units), [], sum(All_signif_units)), 'Centered', true);
fit_both_pca_coeff = pca(reshape(cat(1,fit_FR_init(:,t_pt_for_pca ,All_signif_units),fit_FR_cor(:,t_pt_for_pca ,All_signif_units)), [], sum(All_signif_units)), 'Centered', true);
fit_init_pca_coeff = [fit_init_pca_coeff(:,1:2), fit_both_pca_coeff];
[fit_init_pca_coeff,~] = qr(fit_init_pca_coeff);
fit_init_pca_coeff= fit_init_pca_coeff(:,1:7);
fit_cor_pca_coeff = [fit_cor_pca_coeff(:,1:2), fit_both_pca_coeff];
[fit_cor_pca_coeff,~] = qr(fit_cor_pca_coeff);
fit_cor_pca_coeff= fit_cor_pca_coeff(:,1:7);

score_init_init = permute( pagemtimes( fit_init_pca_coeff', permute(fit_FR_init(:,:,All_signif_units), [3,1,2]) ), [2,3,1]);
score_cor_init = permute( pagemtimes( fit_init_pca_coeff', permute(fit_FR_cor(:,:,All_signif_units), [3,1,2]) ), [2,3,1]);

score_cor_cor = permute( pagemtimes( fit_cor_pca_coeff', permute(fit_FR_cor(:,:,All_signif_units), [3,1,2]) ), [2,3,1]);
score_init_cor = permute( pagemtimes( fit_cor_pca_coeff', permute(fit_FR_init(:,:,All_signif_units), [3,1,2]) ), [2,3,1]);

score_cor_both = permute( pagemtimes( fit_both_pca_coeff', permute(fit_FR_cor(:,:,All_signif_units), [3,1,2]) ), [2,3,1]);
score_init_both = permute( pagemtimes( fit_both_pca_coeff', permute(fit_FR_init(:,:,All_signif_units), [3,1,2]) ), [2,3,1]);

plane_overlap_init = fit_both_pca_coeff'*fit_init_pca_coeff;
plane_overlap_cor = fit_both_pca_coeff'*fit_cor_pca_coeff;
plane_overlap_init(1:3,1) = plane_overlap_init(1:3,1)./norm(plane_overlap_init(1:3,1));
plane_overlap_init(1:3,2) = plane_overlap_init(1:3,2)./norm(plane_overlap_init(1:3,2));
plane_overlap_cor(1:3,1) = plane_overlap_cor(1:3,1)./norm(plane_overlap_cor(1:3,1));
plane_overlap_cor(1:3,2) = plane_overlap_cor(1:3,2)./norm(plane_overlap_cor(1:3,2));

cov_init = cov(reshape(fit_FR_init(:,:,All_signif_units),[],sum(All_signif_units)));
cov_cor = cov(reshape(fit_FR_cor(:,:,All_signif_units),[],sum(All_signif_units)));


overlap = trace(cov_init*cov_cor)./(norm(cov_init,'fro')*norm(cov_cor,'fro'));

tot_init_var = sum(sum(sum( (fit_FR_init(:,:,All_signif_units)-mean(mean(fit_FR_init(:,:,All_signif_units),1),2) ).^2 )));
tot_cor_var = sum(sum(sum( (fit_FR_cor(:,:,All_signif_units)-mean(mean(fit_FR_cor(:,:,All_signif_units),1),2) ).^2 )));
PC12_var_init_init = sum(sum(sum( (score_init_init(:,:,1:2)-mean(mean(score_init_init(:,:,1:2),1),2) ).^2 )));
PC12_var_cor_init = sum(sum(sum( (score_cor_init(:,:,1:2)-mean(mean(score_cor_init(:,:,1:2),1),2) ).^2 )));
PC12_var_init_cor = sum(sum(sum( (score_init_cor(:,:,1:2)-mean(mean(score_init_cor(:,:,1:2),1),2) ).^2 )));
PC12_var_cor_cor = sum(sum(sum( (score_cor_cor(:,:,1:2)-mean(mean(score_cor_cor(:,:,1:2),1),2) ).^2 )));


percent_init_init = 100*PC12_var_init_init/tot_init_var;
percent_cor_init = 100*PC12_var_cor_init/tot_cor_var;
percent_init_cor = 100*PC12_var_init_cor/tot_init_var;
percent_cor_cor = 100*PC12_var_cor_cor/tot_cor_var;

num_signif_units = sum(All_signif_units);
fit_FR_init_signif = fit_FR_init(:,:,All_signif_units);
fit_FR_cor_signif = fit_FR_cor(:,:,All_signif_units);
for k = 1:1000
    curr_units = randi(num_signif_units ,[num_signif_units,1] );
    cov_init = cov(reshape(fit_FR_init_signif(:,:, curr_units),[],num_signif_units ));
    cov_cor = cov(reshape(fit_FR_cor_signif(:,:, curr_units),[],num_signif_units));


    overlap_bs(k) = trace(cov_init*cov_cor)./(norm(cov_init,'fro')*norm(cov_cor,'fro'));

    tot_init_var_bs = sum(sum(sum( (fit_FR_init_signif(:,:,curr_units)-mean(mean(fit_FR_init_signif(:,:,curr_units),1),2) ).^2 )));
    tot_cor_var_bs = sum(sum(sum( (fit_FR_cor_signif(:,:,curr_units)-mean(mean(fit_FR_cor_signif(:,:,curr_units),1),2) ).^2 )));

    fit_init_pca_coeff_bs = pca(reshape(fit_FR_init_signif(:,t_pt_for_pca ,curr_units), [], num_signif_units), 'Centered', true);
    fit_cor_pca_coeff_bs = pca(reshape(fit_FR_cor_signif(:,t_pt_for_pca ,curr_units), [], num_signif_units), 'Centered', true);


    score_init_init_bs = permute( pagemtimes( fit_init_pca_coeff_bs', permute(fit_FR_init_signif(:,:,curr_units), [3,1,2]) ), [2,3,1]);
    score_cor_init_bs = permute( pagemtimes( fit_init_pca_coeff_bs', permute(fit_FR_cor_signif(:,:,curr_units), [3,1,2]) ), [2,3,1]);

    score_cor_cor_bs = permute( pagemtimes( fit_cor_pca_coeff_bs', permute(fit_FR_cor_signif(:,:,curr_units), [3,1,2]) ), [2,3,1]);
    score_init_cor_bs = permute( pagemtimes( fit_cor_pca_coeff_bs', permute(fit_FR_init_signif(:,:,curr_units), [3,1,2]) ), [2,3,1]);




    PC12_var_init_init_bs = sum(sum(sum( (score_init_init_bs(:,:,1:2)-mean(mean(score_init_init_bs(:,:,1:2),1),2) ).^2 )));
    PC12_var_cor_init_bs = sum(sum(sum( (score_cor_init_bs(:,:,1:2)-mean(mean(score_cor_init_bs(:,:,1:2),1),2) ).^2 )));
    PC12_var_init_cor_bs = sum(sum(sum( (score_init_cor_bs(:,:,1:2)-mean(mean(score_init_cor_bs(:,:,1:2),1),2) ).^2 )));
    PC12_var_cor_cor_bs = sum(sum(sum( (score_cor_cor_bs(:,:,1:2)-mean(mean(score_cor_cor_bs(:,:,1:2),1),2) ).^2 )));


    percent_init_init_bs(k) = 100*PC12_var_init_init_bs/tot_init_var_bs;
    percent_cor_init_bs(k) = 100*PC12_var_cor_init_bs/tot_cor_var_bs;
    percent_init_cor_bs(k) = 100*PC12_var_init_cor_bs/tot_init_var_bs;
    percent_cor_cor_bs(k) = 100*PC12_var_cor_cor_bs/tot_cor_var_bs;


    %Calculate plane angle
    fit_both_pca_coeff_bs = pca(reshape(cat(1,fit_FR_init_signif(:,t_pt_for_pca ,curr_units),fit_FR_cor_signif(:,t_pt_for_pca ,curr_units)), [], num_signif_units), 'Centered', true);

    plane_overlap_init_bs = fit_both_pca_coeff_bs'*fit_init_pca_coeff_bs;
    plane_overlap_cor_bs = fit_both_pca_coeff_bs'*fit_cor_pca_coeff_bs;

    init_vec_bs = cross(plane_overlap_init_bs(1:3,1), plane_overlap_init_bs(1:3,2));
    init_vec_bs = init_vec_bs./norm(init_vec_bs);
    cor_vec_bs = cross(plane_overlap_cor_bs(1:3,1), plane_overlap_cor_bs(1:3,2));
    cor_vec_bs = cor_vec_bs./norm(cor_vec_bs);
    plane_angle_bs(k) = acosd(abs(dot( init_vec_bs, cor_vec_bs )));




end

disp(['Overlap: ' num2str(overlap,3), ', ' num2str(prctile(overlap_bs,[2.5,97.5]),3)])
disp(['Percentage Init in Init space: ', num2str(percent_init_init,3), ', ', num2str(prctile(percent_init_init_bs,[2.5,97.5]),3)])
disp(['Percentage Cor in Init space: ', num2str(percent_cor_init,3), ', ', num2str(prctile(percent_cor_init_bs,[2.5,97.5]),3)])
disp(['Percentage Init in Cor space: ', num2str(percent_init_cor,3), ', ', num2str(prctile(percent_init_cor_bs,[2.5,97.5]),3)])
disp(['Percentage Cor in Cor space: ', num2str(percent_cor_cor,3), ', ', num2str(prctile(percent_cor_cor_bs,[2.5,97.5]),3)])


center_init_both = squeeze(mean(score_init_both(:,1,:),1));
center_cor_both = squeeze(mean(score_cor_both(:,1,:),1));
init_plane = center_init_both(1:3) + 200*[plane_overlap_init(1:3,1), plane_overlap_init(1:3,2), -plane_overlap_init(1:3,1), -plane_overlap_init(1:3,2)];
cor_plane = center_cor_both(1:3) + 200*[plane_overlap_cor(1:3,1), plane_overlap_cor(1:3,2), -plane_overlap_cor(1:3,1), -plane_overlap_cor(1:3,2)];
center_init_init = squeeze(mean(score_init_init(:,1,:),1));
init_plane_init = center_init_init(1:3) + 200*[1,1,0; -1,1,0; -1,-1,0; 1,-1,0]';
center_cor_cor = squeeze(mean(score_cor_cor(:,1,:),1));
cor_plane_cor = center_cor_cor(1:3) + 200*[1,1,0; -1,1,0; -1,-1,0; 1,-1,0]';

[tmp,~] = qr(plane_overlap_init(1:3,1:2));
plane_overlap_init(1:3,1:2) = tmp(1:3,1:2);
[tmp,~] = qr(plane_overlap_cor(1:3,1:2));
plane_overlap_cor(1:3,1:2) = tmp(1:3,1:2);

%Calculate plane angle
init_vec = cross(plane_overlap_init(1:3,1), plane_overlap_init(1:3,2));
init_vec = init_vec./norm(init_vec);
cor_vec = cross(plane_overlap_cor(1:3,1), plane_overlap_cor(1:3,2));
cor_vec = cor_vec./norm(cor_vec);
plane_angle = acosd(abs(dot( init_vec, cor_vec )));

disp(['Angle between planes: ', num2str(plane_angle,'%2.1f'), ', ', num2str(prctile(plane_angle_bs,[2.5,97.5]),3)])



%Figure 6, 2D version
% %Trajectory is 300ms until 100ms after peak speed
% figure('PaperPosition', [1,1,11.6, 8.75])
% wysiwyg
% ha = my_subplot([1 1], 1, 3, [.05,.2,.15],.15);
% set(gcf, 'CurrentAxes', ha{1}(1,1))
% plot(score_init_init(:,:,1)', score_init_init(:,:,2)', 'LineWidth', 2.5)
% hold on
% scatter(score_init_init(:,16,1), score_init_init(:,16,2), 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
% axis equal
% axis off
% text(160,-50,[num2str(percent_init_init,3) '%'],'Color', init_color, 'FontSize',18)
% annotation('rectangle', ha{1}(1,1).Position, 'LineWidth', 2.5, 'Color', init_color)
% 
% set(gcf, 'CurrentAxes', ha{1}(1,2))
% plot(score_cor_init(:,:,1)', score_cor_init(:,:,2)', 'LineWidth', 2.5)
% hold on
% scatter(score_cor_init(:,16,1), score_cor_init(:,16,2), 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
% ha{1}(1,2).XLim = ha{1}(1,1).XLim; 
% ha{1}(1,2).YLim = ha{1}(1,1).YLim; 
% axis off
% text(160,-50,[num2str(percent_cor_init,3) '%'],'Color', cor_color, 'FontSize',18)
% annotation('rectangle', ha{1}(1,2).Position, 'LineWidth', 2.5, 'Color', init_color)
% annotation('textbox',[ha{1}(1,2).Position(1), ha{1}(1,2).Position(2)-.04, ha{1}(1,2).Position(3), .04], 'String', 'Initial Subspace, a.u.', 'Color', init_color, 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')
% 
% set(gcf, 'CurrentAxes', ha{1}(1,3))
% plot(score_init_init(:,:,1)', score_init_init(:,:,2)', 'Color', init_color, 'LineWidth', 2.5)
% hold on
% plot(score_cor_init(:,:,1)', score_cor_init(:,:,2)', 'Color', cor_color, 'LineWidth', 1.5)
% scatter(score_init_init(:,16,1), score_init_init(:,16,2), 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
% scatter(score_cor_init(:,16,1), score_cor_init(:,16,2), 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
% ha{1}(1,3).XLim = ha{1}(1,1).XLim; 
% ha{1}(1,3).YLim = ha{1}(1,1).YLim; 
% axis square
% axis equal
% axis off
% annotation('rectangle', ha{1}(1,3).Position, 'LineWidth', 2.5, 'Color', init_color, 'LineWidth', 2.5)
% 
% set(gcf, 'CurrentAxes', ha{1}(2,2))
% plot(score_cor_cor(:,:,1)', score_cor_cor(:,:,2)', 'LineWidth', 2.5)
% hold on
% scatter(score_cor_cor(:,16,1), score_cor_cor(:,16,2), 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
% axis equal
% axis off
% annotation('rectangle', ha{1}(2,2).Position, 'LineWidth', 2.5, 'Color', cor_color)
% annotation('textbox',[ha{1}(2,2).Position(1), ha{1}(2,2).Position(2)-.04, ha{1}(2,2).Position(3), .04], 'String', 'Corrective Subspace, a.u.', 'Color', cor_color, 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')
% text(160,-50,[num2str(percent_cor_cor,3) '%'],'Color', cor_color, 'FontSize',18)
% ha{1}(2,2).XLim = ha{1}(1,1).XLim; 
% ha{1}(2,2).YLim = ha{1}(1,1).YLim; 
% 
% set(gcf, 'CurrentAxes', ha{1}(2,1))
% plot(score_init_cor(:,:,1)', score_init_cor(:,:,2)', 'LineWidth', 2.5)
% hold on
% scatter(score_init_cor(:,16,1), score_init_cor(:,16,2), 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
% ha{1}(2,1).XLim = ha{1}(2,2).XLim; 
% ha{1}(2,1).YLim = ha{1}(2,2).YLim; 
% text(160,-50,[num2str(percent_cor_init,3) '%'],'Color', init_color, 'FontSize',18)
% axis off
% annotation('rectangle', ha{1}(2,1).Position, 'LineWidth', 2.5, 'Color', cor_color)
% 
% set(gcf, 'CurrentAxes', ha{1}(2,3))
% plot(score_init_cor(:,:,1)', score_init_cor(:,:,2)', 'Color', init_color, 'LineWidth', 1.5)
% hold on
% plot(score_cor_cor(:,:,1)', score_cor_cor(:,:,2)', 'Color', cor_color, 'LineWidth', 2.5)
% scatter(score_init_cor(:,16,1), score_init_cor(:,16,2), 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
% scatter(score_cor_cor(:,16,1), score_cor_cor(:,16,2), 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
% ha{1}(2,3).XLim = ha{1}(2,2).XLim; 
% ha{1}(2,3).YLim = ha{1}(2,2).YLim; 
% axis off
% annotation('rectangle', ha{1}(2,3).Position, 'LineWidth', 2.5, 'Color', cor_color)
% % print(gcf, 'InitCorSubspace', '-dpng')






color_order = get(gca,'ColorOrder');
color_order = [color_order; 0.5,0.9, 0.5];

%%Figure 6
%Plot simulated neural trajectories with 3D planes plus 2D projections
%Trajectory is 300ms until 100ms after peak speed
figure('PaperPosition', [1,1,12.6, 8.75])
wysiwyg
ha = my_subplot([1 1], 1, 3, [.15,.05,.05],[.35,.15,.15]);
set(gcf, 'CurrentAxes', ha{1}(1,2))
plot(score_init_init(:,:,1)', score_init_init(:,:,2)', 'LineWidth', 2.5)
set(gca,'ColorOrder', color_order)
hold on
scatter(score_init_init(:,16,1), score_init_init(:,16,2), 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
axis equal
axis off
curr_plane = [ha{1}(1,2).XLim(1)-5, ha{1}(1,2).YLim(1)-5; ...
    ha{1}(1,2).XLim(2)+5, ha{1}(1,2).YLim(1)-5; ...
    ha{1}(1,2).XLim(2)+5, ha{1}(1,2).YLim(2)+5; ...
    ha{1}(1,2).XLim(1)-5, ha{1}(1,2).YLim(2)+5];
patch(curr_plane(:,1), curr_plane(:,2), [-.1,-.1,-.1,-.1], init_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
text(160,-220,[num2str(percent_init_init,3) '%'],'Color', init_color, 'FontSize',18)
% annotation('rectangle', ha{1}(1,2).Position, 'LineWidth', 2.5, 'Color', init_color)

set(gcf, 'CurrentAxes', ha{1}(1,3))
plot(score_cor_init(:,:,1)', score_cor_init(:,:,2)', 'LineWidth', 2.5)
set(gca,'ColorOrder', color_order)
hold on
scatter(score_cor_init(:,16,1), score_cor_init(:,16,2), 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
ha{1}(1,3).XLim = ha{1}(1,2).XLim;
ha{1}(1,3).YLim = ha{1}(1,2).YLim;
axis off
text(160,-220,[num2str(percent_cor_init,3) '%'],'Color', cor_color, 'FontSize',18)
patch(curr_plane(:,1), curr_plane(:,2), [-.1,-.1,-.1,-.1], init_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
% annotation('rectangle', ha{1}(1,3).Position, 'LineWidth', 2.5, 'Color', init_color)
% annotation('textbox',[ha{1}(1,3).Position(1), ha{1}(1,3).Position(2)-.04, ha{1}(1,3).Position(3), .04], 'String', 'Initial Subspace, a.u.', 'Color', init_color, 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')


set(gcf, 'CurrentAxes', ha{1}(2,3))
plot(score_cor_cor(:,:,1)', score_cor_cor(:,:,2)', 'LineWidth', 2.5)
set(gca,'ColorOrder', color_order)
hold on
scatter(score_cor_cor(:,16,1), score_cor_cor(:,16,2), 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color)
axis equal
axis off
% annotation('rectangle', ha{1}(2,3).Position, 'LineWidth', 2.5, 'Color', cor_color)
text(160,-220,[num2str(percent_cor_cor,3) '%'],'Color', cor_color, 'FontSize',18)
ha{1}(2,3).XLim = ha{1}(1,2).XLim;
ha{1}(2,3).YLim = ha{1}(1,2).YLim;
curr_plane = [ha{1}(2,3).XLim(1)-5, ha{1}(2,3).YLim(1)-5; ...
    ha{1}(2,3).XLim(2)+5, ha{1}(2,3).YLim(1)-5; ...
    ha{1}(2,3).XLim(2)+5, ha{1}(2,3).YLim(2)+5; ...
    ha{1}(2,3).XLim(1)-5, ha{1}(2,3).YLim(2)+5];
patch(curr_plane(:,1), curr_plane(:,2), [-.1,-.1,-.1,-.1], cor_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2)

set(gcf, 'CurrentAxes', ha{1}(2,2))
plot(score_init_cor(:,:,1)', score_init_cor(:,:,2)', 'LineWidth', 2.5)
set(gca,'ColorOrder', color_order)
hold on
scatter(score_init_cor(:,16,1), score_init_cor(:,16,2), 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color)
ha{1}(2,2).XLim = ha{1}(2,3).XLim;
ha{1}(2,2).YLim = ha{1}(2,3).YLim;
text(160,-220,[num2str(percent_cor_init,3) '%'],'Color', init_color, 'FontSize',18)
axis off
% annotation('rectangle', ha{1}(2,2).Position, 'LineWidth', 2.5, 'Color', cor_color)
patch(curr_plane(:,1), curr_plane(:,2), [-.1,-.1,-.1,-.1], cor_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2)

annotation('textbox',[ha{1}(1,2).Position(1), ha{1}(1,2).Position(2)+ha{1}(1,2).Position(4)+.02, ha{1}(2,2).Position(3), .04], 'String', 'Initial Submovements', 'Color', init_color, 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')
annotation('textbox',[ha{1}(1,3).Position(1), ha{1}(1,3).Position(2)+ha{1}(1,3).Position(4)+.02, ha{1}(2,3).Position(3), .04], 'String', 'Corrective Submovements', 'Color', cor_color, 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')
annotation('Line',  (ha{1}(1,3).Position(1)-.016)*[1 1], [0.02, 0.92])

set(gcf, 'CurrentAxes', ha{1}(2,1))
axis off

set(gcf, 'CurrentAxes', ha{1}(1,1))
ha{1}(1,1).Position = [0.01, 0.25, 0.35,0.55];
plot3(score_init_both(:,:,1)', score_init_both(:,:,2)', score_init_both(:,:,3)', 'LineWidth', 2.5)
hold on
set(gca,'ColorOrder', color_order)
plot3(score_cor_both(:,:,1)', score_cor_both(:,:,2)', score_cor_both(:,:,3)', 'LineWidth', 2.5)
patch(init_plane(1,:), init_plane(2,:), init_plane(3,:), init_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
patch(cor_plane(1,:), cor_plane(2,:), cor_plane(3,:), cor_color, 'EdgeColor', 'none', 'FaceAlpha', 0.2)
camzoom(1.6)
axis off
annotation('arrow', [0.29,0.39], [0.65, 0.72], 'Color', init_color, 'LineWidth', 3)
annotation('arrow', [0.29,0.39], [0.32,0.25], 'Color', cor_color, 'LineWidth', 3)
annotation('textbox',[0.15, 0.74, ha{1}(2,2).Position(3), .04], 'String', 'Initial Subspace', 'Color', init_color, 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')
annotation('textbox',[0.15, 0.20, ha{1}(2,2).Position(3), .04], 'String', 'Corrective Subspace', 'Color', cor_color, 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')
annotation('textbox',[0.05, 0.35, .04, .04], 'String', [num2str(plane_angle,'%2.1f') char(176)], 'Color', 'k', 'FontSize', 18, 'EdgeColor', 'none', 'HorizontalAlignment','center')

%Draw normal vectors
norm_init = cross(plane_overlap_init(1:3,1), plane_overlap_init(1:3,2));
norm_cor = cross(plane_overlap_cor(1:3,1), plane_overlap_cor(1:3,2));
norm_init = norm_init/norm(norm_init);
norm_cor = norm_cor/norm(norm_cor);

line_origin = [-125, 100, -130];
line_len = 100;
line(line_origin(1)+[0 line_len*norm_init(1)], line_origin(2)+[0  line_len*norm_init(2)], line_origin(3)+[0  line_len*norm_init(3)], 'Color', init_color, 'LineWidth', 2)
line(line_origin(1)+[0 line_len*norm_cor(1)], line_origin(2)+[0  line_len*norm_cor(2)], line_origin(3)+[0  line_len*norm_cor(3)], 'Color', cor_color, 'LineWidth', 2)

%Angle arc
point1 = 0.4*line_len*norm_init';
point2 =  0.4*line_len*norm_cor';
radius = norm(point1);
normal = cross(point1, point2 );
normal = normal / norm(normal);
point2 = cross(normal,point1);
point1 =  point1/norm( point1);
point2 =  point2/norm( point2);
theta = linspace(0, 2 * pi*(plane_angle/360), 100);
circle_points = line_origin  + radius * (cos(theta)' * point1 + sin(theta)' * point2);
plot3(circle_points(:, 1), circle_points(:, 2), circle_points(:, 3), 'Color', 'k', 'LineWidth', 1)

%%Plot initial and corrective reach vectors, (To be moved)
headLength = 15;
headWidth = 15;
hold on
center = [0.12, 0.85];
init_len = 0.06;
for k = 1:size(fit_vectors_init,1)
ah = annotation('arrow', 'Position', [center(1), center(2), init_len*fit_vectors_init(k,4), init_len*fit_vectors_init(k,5)], ...
            'headStyle', 'cback1', ...
            'HeadLength', headLength, ...
            'HeadWidth', headWidth, ...
            'Units', 'centimeters', 'LineWidth', 3, 'Color', init_color);
end
for k = 1:size(fit_vectors_init,1)
ah = annotation('arrow', 'Position', [center(1), center(2), init_len*.8*fit_vectors_init(k,4), init_len*.8*fit_vectors_init(k,5)], ...
            'headStyle', 'none', ...
            'HeadLength', headLength, ...
            'HeadWidth', headWidth, ...
            'Units', 'centimeters', 'LineWidth', 3, 'Color', color_order(k,:));
end


center = [0.12, 0.15];
corr_len = (mean_cor_speed./mean_init_speed)*init_len;

headLength =  (mean_cor_speed./mean_init_speed)*headLength;
headWidth = (mean_cor_speed./mean_init_speed)*headWidth;
for k = 1:size(fit_vectors_init,1)
ah = annotation('arrow', 'Position', [center(1), center(2), corr_len*fit_vectors_init(k,4), corr_len*fit_vectors_init(k,5)], ...
            'headStyle', 'cback1', ...
            'HeadLength', headLength, ...
            'HeadWidth', headWidth, ...
            'Units', 'centimeters', 'LineWidth', 3, 'Color', cor_color);
end
for k = 1:size(fit_vectors_init,1)
ah = annotation('arrow', 'Position', [center(1), center(2), corr_len*.8*fit_vectors_init(k,4), corr_len*.8*fit_vectors_init(k,5)], ...
            'headStyle', 'none', ...
            'HeadLength', headLength, ...
            'HeadWidth', headWidth, ...
            'Units', 'centimeters', 'LineWidth', 3, 'Color', color_order(k,:));
end

% print(gcf, 'Figure6c.svg', '-dsvg', '-vector')
% svg_fix_viewbox('Figure6c.svg', 'Figure6d.svg')

%% print(gcf, 'Figure6.pdf', '-dpdf', '-vector')
%% print(gcf, 'Figure6.png', '-dpng')
%% fig2svg('./figure_graphics/figure6.svg', gcf)





