%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

%% Specify the address of current folder
cd 'D:\~Work\github\d-stem-LUR\code'

%% Clean workspace and load data
clc
clear all
close all
load('..\data\Tehran2015_PM2point5.mat');

%% Remove predictors that cannot let h-block CV computability due to small number of distinct values
feasible_hNum = 5;
predictors = FixPredsPool4hBlockCV(predictors , feasible_hNum);


%% simulation that produce missing, impute them, and assess imputed values
input = response_ground;
n = 50; % specify number of simulation repeats, for example 20
simul_percent = 10 ; % specify percent of added random missing values
sites = 1:length(input.lat);   % sites = 1:30;
% extract data values for subsequent analyses
input.data = input.data(sites,:);
input.site_names = input.site_names(sites,:);
input.nonmiss_count = input.nonmiss_count(sites,:);
input.nonmiss_percent = input.nonmiss_percent(sites,:);
input.lon = input.lon(sites,:);
input.lat = input.lat(sites,:);
% remember orginal data and orginal missing pattern:
input.nonmiss = ~isnan(input.data);
input.nonmiss_freq = sum(input.nonmiss);
input.orginal_data = input.data; % We will add random missing to input.data
% plot orginal missing pattern
[rsize, csize] = size(input.nonmiss); 
imagesc((1:csize)+0.5, (1:rsize)+0.5, input.nonmiss);
colormap(gray);                              % Use a gray colormap
% start of loop
output = nan(n,16);
r = 1;   % change this from 1 to n when you want manually run the involved code of below loop
for r = 1: n
% produce missing:
input.data = input.orginal_data; % restore orginal data in each loop repeat
input.n = sum(input.nonmiss_freq);
input.simul_percent = simul_percent;
input.simul_m = round(input.n*input.simul_percent/100);
count = 0;
nrow = size(input.data,1);
ncol = size(input.data,2);
freq = input.nonmiss_freq;
input.simul_nonmiss = input.nonmiss;
min_freq_rancol = min(freq);
while count < input.simul_m
    rancol = random('unid',ncol);
    if freq(rancol) < min_freq_rancol
        min_freq_rancol = freq(rancol);
    end
    if freq(rancol)>1
        ranrow = random('unid',nrow);
        if input.simul_nonmiss(ranrow,rancol)==1
            input.simul_nonmiss(ranrow,rancol)=0;
            count = count+1;
        end
    end
end
% % plot missing pattern in dataset with added random missing
% [rsize, csize] = size(input.simul_nonmiss); 
% imagesc((1:csize)+0.5, (1:rsize)+0.5, input.simul_nonmiss);
% colormap(gray);                              % Use a gray colormap
% % set(gca, 'XTick', 1:(csize+1), 'YTick', 1:(rsize+1), ...  % Change some axes properties
% %          'XLim', [1 csize+1], 'YLim', [1 rsize+1], ...
% %          'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
input.simulated = (input.simul_nonmiss ~= input.nonmiss); % indicate simulated missing
input.simulated_freq = sum(input.simulated); 
tabulate(input.simulated_freq)% freq of number of induced missing data over time
% create simulated data that contain some random missing generated data:
input.data(input.simulated)=nan;
% impute missing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('..\d-stem\Src')    % path to the d-stem Src folder is required.
index = [62 95 132 58]; % define index of spatial predictors
reffect = []; % define random effect predictors based on thier location in index
cvsites = [];
predicts = predictors;
predicts.X = predictors.X(sites,:);
obj_stem_model_LUR = FitModelSpatialPredictors(input , predicts , index , reffect , cvsites,NaN);
obj_stem_model_krig = FitModelSpatialPredictors(input , predicts , [] , reffect , cvsites,NaN);
obj_stem_model_krig_ACE = FitModelSpatialPredictorsAutoCorrExp(input , predicts , [] , reffect , cvsites,NaN);
% impute by using station average in place of missings
model_meansubstit = repmat(nanmean(input.data), size(input.data,1), 1);
% extract simulated data
sim_observed = input.orginal_data(input.simulated);
sim_model_LUR = obj_stem_model_LUR.stem_EM_result.y_hat_back(input.simulated);
sim_model_krig = obj_stem_model_krig.stem_EM_result.y_hat_back(input.simulated);
sim_model_krig_ACE = obj_stem_model_krig_ACE.stem_EM_result.y_hat_back(input.simulated);
sim_model_meansubstit = model_meansubstit (input.simulated);
% compute power of missing imputation of models
sim_MAE_meansubstit = mean(abs(sim_model_meansubstit-sim_observed)); % mean absolute error (MAE)
sim_MAE_krig_ACE = mean(abs(sim_model_krig_ACE-sim_observed)); % mean absolute error (MAE)
sim_MAE_krig = mean(abs(sim_model_krig-sim_observed)); % mean absolute error (MAE)
sim_MAE_LUR = mean(abs(sim_model_LUR-sim_observed)); % mean absolute error (MAE)
sim_MAPE_meansubstit = 100*mean(abs((sim_model_meansubstit-sim_observed)./sim_observed)); % mean absolute percentage error (MAPE)
sim_MAPE_krig_ACE = 100*mean(abs((sim_model_krig_ACE-sim_observed)./sim_observed)); % mean absolute percentage error (MAPE)
sim_MAPE_krig = 100*mean(abs((sim_model_krig-sim_observed)./sim_observed)); % mean absolute percentage error (MAPE)
sim_MAPE_LUR = 100*mean(abs((sim_model_LUR-sim_observed)./sim_observed)); % mean absolute percentage error (MAPE)
sim_RMSE_meansubstit = sqrt(mean(power(sim_model_meansubstit-sim_observed,2))); % Root Mean Square Error (RMSE)
sim_RMSE_krig_ACE = sqrt(mean(power(sim_model_krig_ACE-sim_observed,2))); % Root Mean Square Error (RMSE)
sim_RMSE_krig = sqrt(mean(power(sim_model_krig-sim_observed,2))); % Root Mean Square Error (RMSE)
sim_RMSE_LUR = sqrt(mean(power(sim_model_LUR-sim_observed,2))); % Root Mean Square Error (RMSE)
sim_PNV_meansubstit = 100*mean(sim_model_meansubstit<0); % Percent of Negative Values (PNV)
sim_PNV_krig_ACE = 100*mean(sim_model_krig_ACE<0); % Percent of Negative Values (PNV)
sim_PNV_krig = 100*mean(sim_model_krig<0); % Percent of Negative Values (PNV)
sim_PNV_LUR = 100*mean(sim_model_LUR<0); % Percent of Negative Values (PNV)
summary = [sim_MAE_meansubstit sim_MAE_krig_ACE sim_MAE_krig sim_MAE_LUR sim_MAPE_meansubstit sim_MAPE_krig_ACE sim_MAPE_krig sim_MAPE_LUR sim_RMSE_meansubstit sim_RMSE_krig_ACE sim_RMSE_krig sim_RMSE_LUR sim_PNV_meansubstit sim_PNV_krig_ACE sim_PNV_krig sim_PNV_LUR ];
output(r,:) = summary;
end

% export the raw results of simulation
varnames ={'MAE_meansubstit', 'MAE_krig_ACE', 'MAE_krig', 'MAE_LUR', 'MAPE_meansubstit', 'MAPE_krig_ACE', 'MAPE_krig', 'MAPE_LUR', 'RMSE_meansubstit', 'RMSE_krig_ACE', 'RMSE_krig', 'RMSE_LUR', 'PNV_meansubstit', 'PNV_krig_ACE', 'PNV_krig', 'PNV_LUR'};
table_output = array2table(round(output,2),'VariableNames',varnames);
writetable(table_output,'output_miss_simul_raw.txt','Delimiter',';');

% export the average of simulation results
outmean = transpose(nanmean(output));
table_outmean = array2table(round(outmean,2),'RowNames',varnames);
writetable(table_outmean,'output_miss_simul_mean.txt','WriteRowNames', true, 'WriteVariableNames', false, 'Delimiter', '\t');
