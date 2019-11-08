%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

function summary = AssessModel( observed_data , modeled_data)
% Assess fitted model using observed_data and modeled_data
% The inputs are n1*T matrixs, where n1 is number of sites and T is length of Time.

% Save observed and modeled data
summary.observed_data = observed_data;
summary.modeled_data = modeled_data;

% Find cells that have nonmissing data and select nonmissed modeled data
summary.nonmissing_data = ~isnan(summary.observed_data);
summary.nonmiss_modeled_data = summary.modeled_data;
summary.nonmiss_modeled_data(~summary.nonmissing_data)= nan;

% Create missing imputed data by replaceing missing by modeled data
summary.missing_imputed = summary.observed_data;
summary.missing_imputed(~summary.nonmissing_data) = summary.modeled_data(~summary.nonmissing_data);

% Report station missing percent and max missing among all stations
summary.station_miss_percent = 1- sum(summary.nonmissing_data,2)/ size(summary.nonmissing_data,2);
summary.station_miss_max = max(summary.station_miss_percent);

% Report total nonmissing count and percent
summary.nonmissing_count = sum(sum(summary.nonmissing_data));
summary.datasetcells_count = size(summary.observed_data, 1) * size(summary.observed_data, 2);
summary.total_miss_percent = 1- (summary.nonmissing_count / summary.datasetcells_count);

% Produce observed and modeled vector data based on nonmissing
summary.observed_nonmiss = summary.observed_data(summary.nonmissing_data);
summary.modeled_nonmiss = summary.modeled_data(summary.nonmissing_data);

% Produce temporal observed and modeled vector data based on nonmissing
summary.tobserved_nonmiss = transpose(nanmean(summary.observed_data));
summary.tmodeled_nonmiss = transpose(nanmean(summary.nonmiss_modeled_data));

% Produce spatial observed and modeled vector data based on nonmissing
summary.sobserved_nonmiss = nanmean(summary.observed_data,2);
summary.smodeled_nonmiss = nanmean(summary.nonmiss_modeled_data,2);

% Calculate criteria
Criteria_name = {'r'; 'MB'; 'MAGE'; 'RMSE'; 'MNB'; 'MNAE'; 'NMB'; 'NMAE'; 'FB'; 'FAE'; 'MNFB'; 'MNAFE'; 'NMBF'; 'NMAEF'; 'NNR'; 'WNNR'; 'CD'; 'R2'; 'MSE'};
summary.Criteria_name = Criteria_name;
Criteria_label = {'Correlation coefficien'; 'Mean Bias'; 'Mean Absolute Gross Error'; 'Root Mean Square Error'; 'Mean Normalized Bias'; 'Mean Normalized Absolute Error'; 'Normalized Mean Bias'; 'Normalized Mean Absolute Error'; 'Fractional Bias'; 'Fractional Absolute Error'; 'Mean Normalized Factor Bias'; 'Mean normalized absolute factor error'; 'Normalized Mean Bias Factor'; 'Normalized Mean Absolute Error Factor'; 'Normalised mean square error of the Normalised Ratios'; 'Weighted Normalised mean square error of the Normalised Ratios'; 'Coefficient of Determination'; 'R-Squared'; 'Mean Square Error'};
summary.Criteria_label = Criteria_label;
summary.Criteria_spatiotemporal = EvaluationCriteria(summary.observed_nonmiss, summary.modeled_nonmiss);
summary.Criteria_spatial = EvaluationCriteria(summary.sobserved_nonmiss, summary.smodeled_nonmiss);
summary.Criteria_temporal = EvaluationCriteria(summary.tobserved_nonmiss, summary.tmodeled_nonmiss);

end
