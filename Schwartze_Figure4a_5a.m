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
%     load([data_location 'Q_Spikes_20180419-data.mat'])
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
neuron_ID{d} = find(M1_units_i);

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

all_x_init{d}=X_init;
all_x_cor{d}=X_cor;

for n = 1:length(best_init_t{d})    
best_SpikeFiringRates_peakVel_init{d}(:,n) = SpikeFiringRates_peakVel_init(:,(20+best_init_t{d}(n)),n);
best_SpikeFiringRates_peakVel_cor{d}(:,n) = SpikeFiringRates_peakVel_cor(:,(20+best_cor_t{d}(n)),n);
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

clearvars -except monkey date_strings d avg_speed_init avg_speed_cor best_R2_tot_init best_init_t best_R2_tot_cor best_cor_t signif_units best_Coefficients_init best_Coefficients_cor neuron_ID ...
    best_std_FR_init best_std_FR_cor pref_dir_init pref_dir_cor best_R2_partial_init best_R2_partial_cor RMSE RMSE_init RMSE_cor best_RMSE_init best_RMSE_cor best_SpikeFiringRates_peakVel_init best_SpikeFiringRates_peakVel_cor all_x_init all_x_cor

d
end
   save(['R:\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results_withfiringrate'])
else
    load(['\\kumc.edu\data\Research\SOM RSCH\RouseLab\DataFiles\Project_Data\20160504_COT_precision\data_analyses\COT_Direction_Regress\' monkey '_regress_results_withfiringrate'])
end

All_day_index = cellfun(@(x) ones(size(x,3),1), signif_units, 'UniformOutput', false);
for d = 1:length(signif_units)
    All_day_index{d} = d*ones(size(signif_units{d},3),1);
end
All_day_index = cat(1,All_day_index{:});
    
    
All_avg_speed_init = avg_speed_init(All_day_index);
All_avg_speed_cor  = avg_speed_cor(All_day_index);

All_signif_units = squeeze(cat(3,signif_units{:}));
% All_SpikeFiringRates_init = cat(2,best_SpikeFiringRates_peakVel_init{:});
% All_SpikeFiringRates_cor = cat(2,best_SpikeFiringRates_peakVel_cor{:});
All_Coefficients_init = cat(2,best_Coefficients_init{:});
All_Coefficients_cor = cat(2,best_Coefficients_cor{:});
All_pref_dir_init = cat(2,pref_dir_init{:});
All_pref_dir_cor  = cat(2,pref_dir_cor{:});
All_best_R2_tot_init = squeeze(cat(3,best_R2_tot_init{:}));
All_best_R2_tot_cor = squeeze(cat(3,best_R2_tot_cor{:}));
All_best_R2_partial_init = cat(2,best_R2_partial_init{:});
All_best_R2_partial_cor  = cat(2,best_R2_partial_cor{:});
All_best_RMSE_init = cat(2,best_RMSE_init{:});
All_best_RMSE_cor = cat(2,best_RMSE_cor{:});

%Velocity regression coefficients depth of modulation
vel_Modulation_init = 2*sqrt(sum(All_Coefficients_init(4:5,:).^2,1));
vel_Modulation_cor  = 2*sqrt(sum(All_Coefficients_cor(4:5,:).^2,1));
    
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


for d = 1:length(best_SpikeFiringRates_peakVel_init)
    Coefficients_init = All_Coefficients_init(:,All_day_index==d);
    Coefficients_cor = All_Coefficients_cor(:,All_day_index==d);
    for n = 1:size(best_SpikeFiringRates_peakVel_init{d},2)
        %Pred spike firing rates
%         pred_SpikeFiringRates_peakVel_init{d}(:,n) =  [ones(size(all_x_init{d},1),1); all_x_init{d}]*Coefficients_init(:,n);
        %Subtract predicted position and speed terms
        marg_SpikeFiringRates_peakVel_init{d}(:,n) =  best_SpikeFiringRates_peakVel_init{d}(:,n) - [all_x_init{d}(:,[1,2,5])]*Coefficients_init([2,3,6],n);
        marg_SpikeFiringRates_peakVel_cor{d}(:,n) =  best_SpikeFiringRates_peakVel_cor{d}(:,n) - [all_x_cor{d}(:,[1,2,5])]*Coefficients_cor([2,3,6],n);
  end
end

for d = 1:length(all_x_init)
    all_vel_dir_init{d} = atan2d(all_x_init{d}(:,4),all_x_init{d}(:,3));
    all_vel_dir_cor{d} = atan2d(all_x_cor{d}(:,4),all_x_cor{d}(:,3));
end


% P_Spikes_20170630-data d=1
% -	Neuron 213
% 
% P_Spikes_20170712-data, d=2
% -	Neuron 64
% 
% P_Spikes_20170703-data, d=3
% -	Neuron 259

angles = -180:1:180;




ex_neuron1 = 178;
    d = All_day_index(ex_neuron1);
    day_index = find(All_day_index==d);
    n = find(day_index==ex_neuron1);
     fr_hat_init = All_Coefficients_init([1,4,5],ex_neuron1)'*[ones(size(angles));avg_speed_init(d)*cosd(angles);avg_speed_init(d)*sind(angles)];
    fr_hat_cor = All_Coefficients_cor([1,4,5],ex_neuron1)'*[ones(size(angles));avg_speed_cor(d)*cosd(angles);avg_speed_cor(d)*sind(angles)];
     
    figure;
scatter(all_vel_dir_init{d}, marg_SpikeFiringRates_peakVel_init{d}(:,n), 10, '+',  'MarkerEdgeColor', init_color+.25)
hold on
scatter(all_vel_dir_cor{d}, marg_SpikeFiringRates_peakVel_cor{d}(:,n), 10, 'x', 'MarkerEdgeColor', cor_color+.25)
legend({'Initial', 'Corrective'})
plot(angles,fr_hat_init, 'Color', init_color,'HandleVisibility','off')
plot(angles,fr_hat_cor, 'Color', cor_color,'HandleVisibility','off')
scatter(pref_dir_init{d}(n), max(fr_hat_init), 30, 'o', 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color,'HandleVisibility','off')
scatter(pref_dir_cor{d}(n), max(fr_hat_cor), 30, 'o', 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color,'HandleVisibility','off')
scatter(wrapTo180(pref_dir_init{d}(n)+180), min(fr_hat_init), 30, 'o', 'MarkerEdgeColor', init_color, 'MarkerFaceColor', [1,1,1],'HandleVisibility','off')
scatter(wrapTo180(pref_dir_cor{d}(n)+180), min(fr_hat_cor), 30, 'o', 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', [1,1,1],'HandleVisibility','off')
xlim([-181,181])
set(gca,'XTick', -180:90:180)
ylabel('Firing Rate, Hz')
xlabel('Angle, degrees')
set(gca, 'FontSize', 18)
% print(gcf, 'ExNeuron1', '-dpdf')
% print(gcf, 'ExNeuron1', '-dtiff')

disp(['Preferred Init Dir: ' num2str(pref_dir_init{d}(n))])
disp(['Preferred Corr Dir: ' num2str(pref_dir_cor{d}(n))])

ex_neuron2 = 785;
    d = All_day_index(ex_neuron2);
    day_index = find(All_day_index==d);
    n = find(day_index==ex_neuron2);
    fr_hat_init = All_Coefficients_init([1,4,5],ex_neuron2)'*[ones(size(angles));avg_speed_init(d)*cosd(angles);avg_speed_init(d)*sind(angles)];
    fr_hat_cor = All_Coefficients_cor([1,4,5],ex_neuron2)'*[ones(size(angles));avg_speed_cor(d)*cosd(angles);avg_speed_cor(d)*sind(angles)];



figure;
scatter(all_vel_dir_init{d}, marg_SpikeFiringRates_peakVel_init{d}(:,n), 10, '+',  'MarkerEdgeColor', init_color+.25)
hold on
scatter(all_vel_dir_cor{d}, marg_SpikeFiringRates_peakVel_cor{d}(:,n), 10, 'x', 'MarkerEdgeColor', cor_color+.25)
legend({'Initial', 'Corrective'})
plot(angles,fr_hat_init, 'Color', init_color,'HandleVisibility','off')
plot(angles,fr_hat_cor, 'Color', cor_color,'HandleVisibility','off')
scatter(pref_dir_init{d}(n), max(fr_hat_init), 30, 'o', 'MarkerEdgeColor', init_color, 'MarkerFaceColor', init_color,'HandleVisibility','off')
scatter(pref_dir_cor{d}(n), max(fr_hat_cor), 30, 'o', 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', cor_color,'HandleVisibility','off')
scatter(wrapTo180(pref_dir_init{d}(n)+180), min(fr_hat_init), 30, 'o', 'MarkerEdgeColor', init_color, 'MarkerFaceColor', [1,1,1],'HandleVisibility','off')
scatter(wrapTo180(pref_dir_cor{d}(n)+180), min(fr_hat_cor), 30, 'o', 'MarkerEdgeColor', cor_color, 'MarkerFaceColor', [1,1,1],'HandleVisibility','off')
xlim([-181,181])
set(gca,'XTick', -180:90:180)
ylabel('Firing Rate, Hz')
xlabel('Angle, degrees')
set(gca, 'FontSize', 18)
% print(gcf, 'ExNeuron2', '-dpdf')
% print(gcf, 'ExNeuron2', '-dtiff')

disp(['Init DOM: ' num2str( max(fr_hat_init)-min(fr_hat_init))])
disp(['Corr DOM: ' num2str(max(fr_hat_cor)-min(fr_hat_cor))])


