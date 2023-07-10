function [R2_partial, R2_tot, Coefficients, p_vals, RMSE] = Do_mutiregress(X,firingrate)
%Input  - X Independent variables, samples x n variable
%       - firing rate for one neuron, samples x 1
%Output - RR_Partial for [position, velocity, speed]


%X=[ones(size(X,1),1) X];
% start_pos = X(:,(1:2));
% vel = X(:,(3:4));
% speed = X(:,5);
%[b,bint,r,rint,stats]=regress(curr_Y(:,1),X);
%[b,~,E,~,~] = mvregress(X,curr_Y);
%(1) intercept (2:3) start_pos (4:5) vel (6) speed]
model=fitlm(X,firingrate);
RRT=1-model.SSR/model.SST;
%ref:https://online.stat.psu.edu/stat462/node/138/
%R2_patial1=SSR(x1,x2,x3)−SSR(x1)/SSE(x1)
%Fit model without Position
model1=fitlm(X(:,3:5),firingrate);
R2_pos = (model.SSR-model1.SSR)/model1.SSE;
%Fit model without Velocity
model2=fitlm([X(:,1:2) X(:,5)],firingrate);
R2_vel = (model.SSR-model2.SSR)/model2.SSE;
%Fit model without Speed
model3=fitlm(X(:,1:4),firingrate);
R2_speed = (model.SSR-model3.SSR)/model3.SSE;

Coefficients = model.Coefficients{:,1}; %(1) intercept (2:3) start_pos (4:5) vel (6) speed]
p_vals = model.Coefficients{:,4};

rho=partialcorri(firingrate,X);%[start_pos vel speed]
R2_partial=rho.^2;
R2_partial = [sum(R2_partial(1:2)), sum(R2_partial(3:4)), sum(R2_partial(5))];
R2_tot = model.Rsquared.ordinary;
RMSE = model.RMSE;
end


