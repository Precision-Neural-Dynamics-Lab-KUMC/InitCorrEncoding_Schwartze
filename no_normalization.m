function [X, SD] = no_normalization(start_pos, end_pos, timestep)
%X = normalization(start_pos, end_pos, timestep)
%Calculates (starting) position, velocity, and speed of a submovement
%Inputs: start_pos - starting position, 2 columns, x/y
%        end_pos - starting position, 2 columns, x/y
%        timestep - time between start and end, in seconds
%Outputs: X  - values, 5 columns: start_pos (x,y), velocity (dx, dy), speed
%         SD - std. dev., 3 coulumns, start_pos, velocity, speed


vel=(end_pos-start_pos)/timestep;%sec
speed=sqrt(sum(vel.^2,2));
%norm (testData - trainingMean) / trainingStdDev
%X_avg = zeros(1,5);
%pos_SD = std(start_pos);
% X = [start_pos vel speed];

% vel_SD = std(reshape(vel(:,(1:2)),[],1));
% pos_SD = std(reshape(start_pos(:,(1:2)),[],1));
% speed_SD = std(speed(:));
% vel_SD = sqrt(sum(reshape(vel(:,(1:2)),[],1).^2)./numel(vel(:,(1:2))));  %Don't want to substract mean
% pos_SD = sqrt(sum(reshape(start_pos(:,(1:2)),[],1).^2)./numel(start_pos(:,(1:2))));
% speed_SD = sqrt(sum(speed(:).^2)./numel(speed));
vel_SD = mean(sqrt(sum(vel(:,1:2).^2,2)));
pos_SD = mean(sqrt(sum(start_pos(:,1:2).^2,2)));
speed_SD = sqrt(sum(speed(:).^2)./numel(speed));


X = [start_pos vel speed];
SD = [pos_SD, vel_SD, speed_SD];
end
