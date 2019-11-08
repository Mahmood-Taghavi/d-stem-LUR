%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

function predUpdate = fixPredsPool4hBlockCV(predictors , hNum)
% remove predictors which are not suitable for calculation of hBlookCV
% hNum=0 can be used for Leave-One-Out Cross-Validation (LOOCV)

acceptCount = 0;
%i = 1;
for i = 1:size(predictors.X,2)
freq = tabulate(predictors.X(:,i));
freqsize = size(freq,1);
if freqsize > (1+hNum+1)
    acceptCount = acceptCount + 1; 
    predUpdate.sign(acceptCount) = predictors.sign(i);
    predUpdate.sign_names(acceptCount) = predictors.sign_names(i);
    predUpdate.sign_labels(acceptCount) = predictors.sign_labels(i);
    predUpdate.X(: , acceptCount) = predictors.X(: , i);
    predUpdate.X_names(acceptCount) = predictors.X_names(i);
    predUpdate.match(acceptCount) = predictors.match(i);
else
    display(strcat('removedIndex: ',int2str(i),' ; removedName: ',predictors.sign_names(i),' ; removedLabel: ',predictors.sign_labels(i)))
end % end if
end % end for

end % end function
