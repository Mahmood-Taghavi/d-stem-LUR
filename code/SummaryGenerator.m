%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

function summary = SummaryGenerator( obj_stem_model)
% Assess D-Stem spatio-temporal fitted model
% calculations are based on two below articles: 
% Yu, S., et al. (2006). "New unbiased symmetric metrics for evaluation of air quality models." Atmospheric Science Letters 7(1): 26-34.
% Poli, A. A. and M. C. Cirillo (1993). "On the use of the normalized mean square error in evaluating dispersion model performance." Atmospheric Environment. Part A. General Topics 27(15): 2427-2434.
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

% Retrive observed and modeled data
observed_matrix = obj_stem_model.stem_EM_result.y_back;
modeled_matrix = obj_stem_model.stem_EM_result.y_hat_back;

summary = AssessModel(observed_matrix , modeled_matrix);

% cross-validation: Retrive observed and modeled data
if obj_stem_model.cross_validation
    observed_matrix = obj_stem_model.stem_data.stem_crossval.stem_crossval_result{1,1}.y_back;
    modeled_matrix = obj_stem_model.stem_data.stem_crossval.stem_crossval_result{1,1}.y_hat_back;
    summary.cross_validation = AssessModel(observed_matrix , modeled_matrix);
end

end
