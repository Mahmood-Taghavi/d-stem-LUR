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

%% Fit model
addpath('..\d-stem\Src')    % path to the d-stem Src folder is required.

index = []; % [] is empty which means no predictors be used in the kriging model
reffect = []; % define random effect predictors based on thier location in index
cvsites = []; % index of cross validation sites
krig_stem_model = FitModelSpatialPredictors(response_ground , predictors , index , reffect , cvsites,NaN);
krig_summary = SummaryGenerator(krig_stem_model); % See Criteria_spatial and so on inside of krig_summary

index = [62 95 132 58]; % [62 95 132 58] index of spatial predictors in our LUR model
lur_stem_model = FitModelSpatialPredictors(response_ground , predictors , index , reffect , cvsites,NaN);
predictors.X_names(index)
predictors.sign_labels(index)
lur_summary = SummaryGenerator(lur_stem_model); % See Criteria_spatial and so on inside of lur_summary

% Note: the smoothing parameter of Matern autocorrelation function can be set
% in FitModelSpatialPredictors function by choosing between: 'exponential'
% 'matern32' 'matern52'. Currently we set it to 'matern52' (best in our case).

%% Calculate h-block cross-validation metrics
hNum = 2; % define block size in h-block cross-validation, for example 2 or 5, note hNum = 0 lead to LOOCV

index = []; % [] is empty which means no predictors be used in the kriging model
krig_summary_hBlockCV = HBlockCvSpatialPredictors (response_ground , predictors , index , reffect , cvsites, NaN, hNum);
% See hBlockCV inside of krig_summary_hBlockCV, especially the Criteria_spatial and so on inside of it

index = [62 95 132 58]; % [62 95 132 58] index of spatial predictors in our LUR model
lur_summary_hBlockCV = HBlockCvSpatialPredictors (response_ground , predictors , index , reffect , cvsites, NaN, hNum);
% See hBlockCV inside of lur_summary_hBlockCV, especially the Criteria_spatial and so on inside of it
