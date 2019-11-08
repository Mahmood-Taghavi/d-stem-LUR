%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

function output = EvaluationCriteria( observed_data , modeled_data )
% assess D-Stem spatio-temporal fitted model

%name = {'r'; 'MB'; 'MAGE'; 'RMSE'; 'MNB'; 'MNAE'; 'NMB'; 'NMAE'; 'FB'; 'FAE'; 'MNFB'; 'MNAFE'; 'NMBF'; 'NMAEF'; 'NNR'; 'WNNR'; 'CD'; 'R2'; 'MSE'};
%label = {'Correlation coefficien'; 'Mean Bias'; 'Mean Absolute Gross Error'; 'Root Mean Square Error'; 'Mean Normalized Bias'; 'Mean Normalized Absolute Error'; 'Normalized Mean Bias'; 'Normalized Mean Absolute Error'; 'Fractional Bias'; 'Fractional Absolute Error'; 'Mean Normalized Factor Bias'; 'Mean normalized absolute factor error'; 'Normalized Mean Bias Factor'; 'Normalized Mean Absolute Error Factor'; 'Normalised mean square error of the Normalised Ratios'; 'Weighted Normalised mean square error of the Normalised Ratios'; 'Coefficient of Determination'; 'R-Squared'; 'Mean Square Error'};

output = nan(18,1);

% calculate correlation between modeled and observed data
output(1) = corr(modeled_data , observed_data ,'rows','pairwise');

% calculate Mean Bias (MB), mean absolute gross error (MAGE), root mean square error (RMSE)
difference_nonmiss = modeled_data - observed_data;
nonmissing_data = ~isnan(observed_data);
nonmissing_count = sum(sum(nonmissing_data));
output(2) = sum(difference_nonmiss) / nonmissing_count;
output(3) = sum(abs(difference_nonmiss)) / nonmissing_count;
output(4) = sqrt(sum(power(difference_nonmiss,2)) / nonmissing_count);

% calculate Mean Normalized Bias (MNB), Mean normalized absolute error (MNAE)
output(5) = sum(difference_nonmiss ./ observed_data) / nonmissing_count;
output(6) = sum(abs(difference_nonmiss) ./ observed_data) / nonmissing_count;

% calculate Normalized mean bias (MNB), Normalized mean absolute error (NMAE)
m_bar = sum (modeled_data) / nonmissing_count;
o_bar = sum (observed_data) / nonmissing_count;
output(7) = (m_bar / o_bar) -1;
output(8) = (output(3) / o_bar);

% calculate Fractional bias (FB), Fractional absolute error (FAE)
sum_observed_modeled = modeled_data + observed_data;
output(9) = sum(difference_nonmiss ./ sum_observed_modeled) / (2*nonmissing_count);
output(10) = sum(abs(difference_nonmiss) ./ sum_observed_modeled) / (2*nonmissing_count);

% Calculate mean normalized factor bias (MNFB), mean normalized absolute factor error (MNAFE)
g1 = (modeled_data ./ observed_data) -1;
g2 = 1- (observed_data ./ modeled_data);
check = (modeled_data < observed_data);
g = g1;
g(check) = g2(check);
output(11) = sum(g)/nonmissing_count;
output(12) = sum(abs(g))/nonmissing_count;

% Calculate Normalized Mean Bias Factor (NMBF), Normalized mean absolute error factor (NMAEF)
output(13) = output(7);
output(14) = output(3) / o_bar;
if m_bar < o_bar
    output(13) = 1-(o_bar / m_bar); 
    output(14) = output(3) / m_bar;
end

% Calculate (Weighted) Normalised mean square error of the Normalised Ratios
weight = observed_data / o_bar;
normalised_ratio = exp( -1 * abs( log (modeled_data ./ observed_data) ) ); % Elementwise operations using ./ instead of matrix /
output(15) = sum(power(1-normalised_ratio,2)) / sum(normalised_ratio);
output(16) = sum(power(weight,2) .* power(1-normalised_ratio,2)) / sum(weight .* normalised_ratio); % Elementwise operations using .*

% calculate Coefficient of determination
%SSR = sum( (modeled_data - o_bar) .^2 );
SSE = sum( (observed_data - modeled_data) .^2 );
SSTO = sum( (observed_data - o_bar) .^2 );
output(17) = 1- SSE/SSTO;

% calculate Coefficient of determination R-Squared
output(18) = sign(output(1))*output(1)^2;

% calculate Mean Square Error (MSE)
output(19) = output(4)^2;


end
