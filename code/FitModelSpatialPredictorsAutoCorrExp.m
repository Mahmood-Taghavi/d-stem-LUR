%% Authorship and License
% This file is part of the d-stem-LUR project on gitHub.
% https://github.com/Mahmood-Taghavi/d-stem-LUR
% License: GPL v2; Author: Seyed Mahmood Taghavi-Shahri
% Please cite the d-stem-LUR paper in addition to d-stem paper in your work:
% Taghavi-Shahri SM, Fasso A, Mahaki B, Amini H. Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space Time Expectation Maximization.

function obj_stem_model = FitModelSpatialPredictorsAutoCorrExp(response_ground , predictors , index , reffect , cvsites, initial)
% Fit stem model using response ground, and predictors with mentioned index.
% Note: if index=[] then a spatiotemporal kriging model will be fit, 
% reffect specify random effect predictors based on thier location in index
% reffect=[] fit a model with just an intercept random effect term 
% Disable crossvalidation if cvsites set to [] or enable it by using cvsites
% Note: cvsites is index of cross validation sites.
% The initial can be Nan, or an initial vector, or a obj_stem_model.stem_par
% The format of initial vector is [theta_p alpha_p G sigma_eta sigma_eps beta]
% The initial vector can be just beginning elements of the above mentioned vecror
% Note: default autocorrelation function is set in line of "obj_stem_par = stem_par(obj_stem_data,'exponential'); % options: 'exponential' 'matern32' 'matern52'    "

v = ver;
test_stat=any(strcmp('Statistics Toolbox', {v.Name}));
if not(test_stat)
   test_stat=any(strcmp('Statistics and Machine Learning Toolbox', {v.Name}));
end
if not(test_stat)
   error('You need the Statistics Toolbox to run this case study. The Statistics Toolbox seems to be missing.');
end

test_opt=any(strcmp('Optimization Toolbox', {v.Name}));
if not(test_opt)
    error('You need the Optimization Toolbox to run this case study. The Optimization Toolbox seems to be missing.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Data  building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the ground level observations
ground.Y{1} = response_ground.data;
%ground.Y_name{1} = 'PM 2.5';
 ground.Y_name{1} = response_ground.var_names; % auto assign name!
n1 = size(ground.Y{1}, 1);
T = size(ground.Y{1}, 2);

%X_beta
xones = ones(n1,1);
xones_name = {'constant'};
%%with intercept:
ground.X_beta{1} = [xones predictors.X(: , index)];
ground.X_beta_name{1} = [xones_name predictors.X_names(: , index)];
%%without intercept:
% ground.X_beta{1} = predictors.X(: , index);
% ground.X_beta_name{1} = predictors.X_names(: , index);


ground.X_z{1} = ones(n1,1);
% ground.X_z{1} = zeta;
ground.X_z_name{1} = {'constant'};

temp = [xones predictors.X(: , index(reffect))];
temp=reshape(temp,[size(temp,1),1,1,size(temp,2)]);
ground.X_p{1} = temp;
ground.X_p_name{1} = [xones_name predictors.X_names(: , index(reffect))];


% if sum(size(reffect)) == 0
%     obj_stem_varset_p = stem_varset(ground.Y,ground.Y_name,[],[],ground.X_beta,ground.X_beta_name,ground.X_z,ground.X_z_name,[],[]);   
% else
    obj_stem_varset_p = stem_varset(ground.Y,ground.Y_name,[],[],ground.X_beta,ground.X_beta_name,ground.X_z,ground.X_z_name,ground.X_p,ground.X_p_name);   
% end
% 
% if sum(size(index)) == 0
%     obj_stem_varset_p = stem_varset(ground.Y,ground.Y_name,[],[],[],[],ground.X_z,ground.X_z_name,[],[]);   
% end

%coordinates
obj_stem_gridlist_p = stem_gridlist();
ground.coordinates{1} = [response_ground.lat, response_ground.lon];
obj_stem_grid = stem_grid(ground.coordinates{1},'deg','sparse','point');
obj_stem_gridlist_p.add(obj_stem_grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Model building     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj_stem_datestamp = stem_datestamp('01-01-2015 00:00','31-12-2015 00:00',T);

%stem_data object creation
% if test_map
%     shape = shaperead('../Maps/Tehran_Boundary');
% else
    shape=[];
% end

% cross validation using index of sites
if sum(size(cvsites)) == 0
    obj_stem_data = stem_data(obj_stem_varset_p,obj_stem_gridlist_p,[],[],obj_stem_datestamp,shape);
else
    %obj_stem_crossval=stem_crossval({'PM 2.5'},{'point'},{cvsites});
    obj_stem_crossval=stem_crossval({response_ground.var_names},{'point'},{cvsites});
    obj_stem_data = stem_data(obj_stem_varset_p,obj_stem_gridlist_p,[],[],obj_stem_datestamp,shape,[],obj_stem_crossval);    
end

%stem_par object creation
if isa(initial, 'stem_par')
    obj_stem_par = initial;
else
    obj_stem_par = stem_par(obj_stem_data,'exponential'); % options: 'exponential' 'matern32' 'matern52'    
end

%stem_model object creation
obj_stem_model = stem_model(obj_stem_data,obj_stem_par);


%Data transform
obj_stem_model.stem_data.log_transform; % comment to do not use log transform
obj_stem_model.stem_data.standardize;


%obj_stem_par object initialization
if ~isa(initial, 'stem_par')
                            %initial: D-STEM, Test, exp, mat32, mat52
    theta = 12.95;          %such as: 31.42, 100, 21.26, 15.75, 12.95, ...
    alpha = 0.368;          %such as: 0.590, 1, 0.405, 0.390, 0.368, ...
    G = 0.76;               %such as: 0.95, 0.6, 0.77, 0.77, 0.76, ...
    sigma_eta = 0.23;       %such as: 0.03, 0.4, 0.22, 0.22, 0.23, ...
    sigma_eps = 0.295;      %such as: 0.392, 0.1, 0.282, 0.294, 0.295, ...
    
    if and(isa(initial, 'double'),~all(isnan(initial)))
        if length(initial)>=1
            if ~isnan(initial(1))
                theta = initial(1);
            end
        end    
        if length(initial)>=2
            if ~isnan(initial(2))
                alpha = initial(2);
            end
        end           
        if length(initial)>=3
            if ~isnan(initial(3))
                G = initial(3);
            end
        end           
        if length(initial)>=4
            if ~isnan(initial(4))
                sigma_eta = initial(4);
            end
        end           
        if length(initial)>=5
            if ~isnan(initial(5))
                sigma_eps = initial(5);
            end
        end      
    end
    
    obj_stem_par.theta_p = repelem(theta, length(reffect)+1); %km
    obj_stem_par.alpha_p = repelem(alpha, length(reffect)+1); %equall res
    obj_stem_par.G = G;
    obj_stem_par.sigma_eta = sigma_eta;
    obj_stem_par.sigma_eps = sigma_eps;     
    
    %obj_stem_par.v_p = 1
    a = repelem(1, length(reffect)+1);
    b=reshape(a,[1,1,size(a,2)]);
    obj_stem_par.v_p = b;
    
    if and(isa(initial, 'double'),length(initial)==5+length(index)+1)
        if ~any(isnan(initial(6:end)))
           obj_stem_par.beta = initial(6:end);
        else
            obj_stem_par.beta = obj_stem_model.get_beta0();
            userbeta = initial(6:end);
            obj_stem_par.beta(~isnan(userbeta))= userbeta(~isnan(userbeta));
        end
    else
        obj_stem_par.beta = obj_stem_model.get_beta0();
    end
    
%     % custom initial values for beta fixed effect coef:
%     obj_stem_par.beta = 0.019;   % beta initial vector
    
end

obj_stem_model.set_initial_values(obj_stem_par);

%Model estimation
exit_toll = .001;      %default: 0.001, others: 0.005, 0.00001, ...
max_iterations = 100;    %default: 100, others: 50, 500, ...
obj_stem_EM_options = stem_EM_options(exit_toll,max_iterations);
obj_stem_model.EM_estimate(obj_stem_EM_options);
obj_stem_model.set_varcov;  % disable this line to speed-up
obj_stem_model.set_logL;    % disable this line to speed-up

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See results of this model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj_stem_model.print;    
%obj_stem_model.stem_EM_result.stem_kalmansmoother_result.plot;

end
