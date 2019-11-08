%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

function accept_sign = CheckSign( obj_stem_model , predictors , index)
% Check sign of estimated coefficents
if (length(index)) > 0
beta = obj_stem_model.stem_EM_result.stem_par.beta();
beta = beta(2:end);   % remove intercept! Check model involve intercept or not!
beta = transpose(beta); % transpose to agree dimention with expected_sign
if (size(beta,2)~=length(index))
    error('Dimention must agree. Check model involve intercept or not.');
end
expected_sign = predictors.sign(index);
% set arbbitary signs (zero's) as sign beta
expected_sign(expected_sign==0) = sign(beta(expected_sign==0));
accept_sign = ( sum(sign(beta)==expected_sign) == length(expected_sign) );
else
    accept_sign = 1; 
end
