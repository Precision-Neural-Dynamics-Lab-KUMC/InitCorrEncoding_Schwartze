%% 
% Load data

clear variables
% monkey = 'P';
% date_strings = {'20170630', '20170712', '20170703', '20170713', '20170720', '20170731', '20170705', '20170706', '20170714', '20170717', '20170801', '20170802'}; 

monkey = 'Q';
date_strings = {'20180425', '20180426', '20180509', '20180510', '20180529', '20180530', '20180418', '20180419', '20180503', '20180507', '20180619', '20180620'};

% Regression and data organization
run_analysis = true;


if run_analysis

    numCores = feature('numcores');
pool = parpool(numCores-1);

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
    if ~exist (data_location)
        data_location  = '\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_extracted\monk_q\COT_SpikesCombined\';
    end
    load([data_location 'Q_Spikes_20180419-data.mat'])
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
    parfor t = 1:41
    [R2_partial(:,t,n), R2_tot(:,t,n), Coefficients(:,t,n), p_vals(:,t,n), curr_RMSE(t,n)] = Do_mutiregress(X,SpikeFiringRates_peakVel(:,(20+t),n));
    end
end
disp(['Day ' num2str(d) ' all regression done'])
% 
%initial peaks only
num_m = 1000;
R2_partial_init = NaN(3,41,size(SpikeFiringRates_peakVel,3));
R2_tot_init = NaN(1,41,size(SpikeFiringRates_peakVel,3));
Coefficients_init = NaN(6,41,size(SpikeFiringRates_peakVel,3));
Coefficients_init1 = NaN(6,1,size(SpikeFiringRates_peakVel,3),num_m);
Coefficients_init2 = NaN(6,1,size(SpikeFiringRates_peakVel,3),num_m);
p_vals_init = NaN(6,41,size(SpikeFiringRates_peakVel,3));
for n = 1:size(SpikeFiringRates_peakVel,3)
    parfor t = 1:41
    [R2_partial_init(:,t,n), R2_tot_init(:,t,n), Coefficients_init(:,t,n), p_vals_init(:,t,n), curr_RMSE_init(t,n)] = Do_mutiregress(X_init,SpikeFiringRates_peakVel_init(:,(20+t),n));
    end
    [~, best_curr_t] = max(R2_tot_init(:,:,n),[],2);
    parfor m=1:num_m
        rand_trials = randperm(size(X_init,1),size(X_init,1));
        rand_trials1 = rand_trials(1:floor(length(rand_trials)./2));
        rand_trials2 = rand_trials((floor(length(rand_trials)./2)+1):end);
        [~, ~, Coefficients_init1(:,1,n,m), ~, ~] = Do_mutiregress(X_init(rand_trials1,:),SpikeFiringRates_peakVel_init(rand_trials1,(20+best_curr_t),n));
        [~, ~, Coefficients_init2(:,1,n,m), ~, ~] = Do_mutiregress(X_init(rand_trials2,:),SpikeFiringRates_peakVel_init(rand_trials2,(20+best_curr_t),n));
    end
end
disp(['Day ' num2str(d) ' init regression done'])
% 
%corrective peaks only
num_m = 1000;
R2_partial_cor = NaN(3,41,size(SpikeFiringRates_peakVel,3));
R2_tot_cor = NaN(1,41,size(SpikeFiringRates_peakVel,3));
Coefficients_cor = NaN(6,41,size(SpikeFiringRates_peakVel,3));
Coefficients_cor1 = NaN(6,1,size(SpikeFiringRates_peakVel,3),num_m);
Coefficients_cor2 = NaN(6,1,size(SpikeFiringRates_peakVel,3),num_m);
p_vals_cor = NaN(6,41,size(SpikeFiringRates_peakVel,3));
for n = 1:size(SpikeFiringRates_peakVel,3)
    parfor t = 1:41
    [R2_partial_cor(:,t,n), R2_tot_cor(:,t,n), Coefficients_cor(:,t,n), p_vals_cor(:,t,n), curr_RMSE_cor(t,n)] = Do_mutiregress(X_cor,SpikeFiringRates_peakVel_cor(:,(20+t),n));
    end
    [~, best_curr_t] = max(R2_tot_cor(:,:,n),[],2);
    parfor m=1:num_m
        rand_trials = randperm(size(X_cor,1),size(X_cor,1));
        rand_trials1 = rand_trials(1:floor(length(rand_trials)./2));
        rand_trials2 = rand_trials((floor(length(rand_trials)./2)+1):end);
        [~, ~, Coefficients_cor1(:,1,n,m), ~, ~] = Do_mutiregress(X_cor(rand_trials1,:),SpikeFiringRates_peakVel_cor(rand_trials1,(20+best_curr_t),n));
        [~, ~, Coefficients_cor2(:,1,n,m), ~, ~] = Do_mutiregress(X_cor(rand_trials2,:),SpikeFiringRates_peakVel_cor(rand_trials2,(20+best_curr_t),n));
    end
end
disp(['Day ' num2str(d) ' cor regression done'])

RMSE{d} = curr_RMSE;
RMSE_init{d} =  curr_RMSE_init;
RMSE_cor{d} = curr_RMSE_cor;
% 
std_FR_init = NaN(41,size(SpikeFiringRates_peakVel,3));
 std_FR_cor = NaN(41,size(SpikeFiringRates_peakVel,3));
for n = 1:size(SpikeFiringRates_peakVel,3)
    parfor t = 1:41
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
best_p_vals_init{d}(:,n) = p_vals_init(:,best_init_t{d}(n),n);
best_p_vals_cor{d}(:,n) = p_vals_cor(:,best_cor_t{d}(n),n);

end
best_Coefficients_init1{d} = Coefficients_init1;
best_Coefficients_init2{d} = Coefficients_init2;
pref_dir_init1{d} = atan2d(best_Coefficients_init1{d}(5,1,:,:), best_Coefficients_init1{d}(4,1,:,:));
pref_dir_init2{d} = atan2d(best_Coefficients_init2{d}(5,1,:,:), best_Coefficients_init2{d}(4,1,:,:));

best_Coefficients_cor1{d} = Coefficients_cor1;
best_Coefficients_cor2{d} = Coefficients_cor2;
pref_dir_cor1{d} = atan2d(best_Coefficients_cor1{d}(5,1,:,:), best_Coefficients_cor1{d}(4,1,:,:));
pref_dir_cor2{d} = atan2d(best_Coefficients_cor2{d}(5,1,:,:), best_Coefficients_cor2{d}(4,1,:,:));

clearvars -except monkey date_strings d avg_speed_init avg_speed_cor best_R2_tot_init best_init_t best_R2_tot_cor best_cor_t signif_units best_Coefficients_init best_Coefficients_cor ...
    best_std_FR_init best_std_FR_cor pref_dir_init pref_dir_cor best_R2_partial_init best_R2_partial_cor RMSE RMSE_init RMSE_cor best_RMSE_init best_RMSE_cor best_Coefficients_init1 best_Coefficients_init2 best_Coefficients_cor1 best_Coefficients_cor2 best_p_vals_init best_p_vals_cor pref_dir_init1 pref_dir_init2 pref_dir_cor1 pref_dir_cor2
    

d
end
   save(['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results_StatsFinal'])
else
    load(['\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results_StatsFinal'])
end
