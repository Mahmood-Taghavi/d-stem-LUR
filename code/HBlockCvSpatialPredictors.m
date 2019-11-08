%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

function output = hBlockCVdefault (response_ground , predictors , index , reffect , cvsites, initial, hNum)
%%%function output = hBlockCV (response_ground , predictors , index , reffect , cvsites, weather, windex, hNum)
% h-block cross-validation for spatial dependent data to assess fit of stem model 
% input: response ground, pool of predictors, index in pool of predictors, ...
% Note: if index=[] then a spatiotemporal kriging model will be fit, 
% reffect specify random effect predictors based on thier location in index
% reffect=[] fit a model with just an intercept random effect term 
% Disable crossvalidation if cvsites set to [] or enable it by using cvsites
% Note: cvsites is index of cross validation sites.
% hNum is the number of nearest sites that must take apart beside each site
% in the h-block cross-validation. hNum = 0 is equivalent to LOOCV method. 

hNum = round(hNum);
if hNum > length(response_ground.lon)
   hNum = 0; 
end
output.hNum = hNum;

% calculate the distance between sites: 
% dist = obj_stem_model.stem_data.DistMat_p;
lon = response_ground.lon;
lat = response_ground.lat;
dist = nan(length(lon));
for a=1:length(lon)
    xa=lat(a);
    ya=lon(a);
    for b=1:length(lon)
        xb=lat(b);
        yb=lon(b);
        d = sqrt((xa-xb)^2+(ya-yb)^2);
        dist(a,b) = d;
    end
end
output.dist = dist;
output.Basic.response_ground = response_ground;
output.Basic.predictors = predictors;
output.Basic.index = index;
output.Basic.reffect = reffect;
output.Basic.cvsites = cvsites;

output.Basic.obj_stem_model = FitModelSpatialPredictors(response_ground , predictors , index , reffect , cvsites,initial);
output.Basic.accept_sign = CheckSign( output.Basic.obj_stem_model , predictors , index);
output.Basic.summary = SummaryGenerator( output.Basic.obj_stem_model); %!!

observed_matrix = nan(size(output.Basic.obj_stem_model.stem_EM_result.y_back));
modeled_matrix = nan(size(output.Basic.obj_stem_model.stem_EM_result.y_back));
i=1;
repeat = 1;
accept_sign = 1;
for c = 1:size(predictors.X,1)
    if sum(cvsites==c) > 0
        % jump becuase this is one of cvsites
    else
        % find h-nearest sites to this site
        h_nearest = [];
        if hNum>0
            distc = dist(c,:);
            distc(c) = Inf;
            distc(cvsites) = Inf;
            h_nearest = nan(1,hNum);
            for k=1:hNum
                [min_dist,min_index] = min(distc);
                h_nearest(k) = min_index;
                distc(min_index) = Inf;
            end
        end
        h_nearest
        % add c and h_nearest to cvsites for doing h-block cross-validation
        cvsitesplus = [c h_nearest cvsites];
        obj_stem_model = FitModelSpatialPredictors(response_ground , predictors , index , reffect , cvsitesplus,initial);
		accept_sign = accept_sign * CheckSign( obj_stem_model , predictors , index);
        modeled_matrix(i,:)=obj_stem_model.stem_data.stem_crossval.stem_crossval_result{1,1}.y_hat_back(1,:);
        observed_matrix(i,:)=obj_stem_model.stem_data.stem_crossval.stem_crossval_result{1,1}.y_back(1,:);
        i=i+1;
        repeat = repeat + 1;
    end   
end

output.observed_matrix = observed_matrix;
output.modeled_matrix = modeled_matrix;

output.hBlockCV = AssessModel(observed_matrix , modeled_matrix);
output.accept_sign = accept_sign;

end
