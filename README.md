# d-stem-LUR
This is the companion repository for the paper titled "Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space-Time Expectation Maximization" which published in Atmospheric Environment journal (http://doi.org/10.1016/j.atmosenv.2019.117202). 
Currently, we have added the dataset and cleaned code to fit and evaluate the final D-STEM kriging and D-STEM LUR models. Also, the code used for assessment of the imputation models is added.
## data
This folder contains daily PM2.5 concentrations in 2015 and the corresponding pool of 210 potentially predictive variables (PPVs) in 30 monitoring stations in Tehran.
## d-stem
This folder contains the original d-stem software source code, especially its "Src" folder which is required for model building.
## code
This folder contains scripts and functions we are coded to fit D-STEM spatiotemporal kriging and LUR models and also it contains the code developed to assess D-STEM fitted models using various metrics (19 metrics) and from different aspects (spatiotemporal, spatial, and temporal) in h-block cross-validation. Furthermore, the developed code for assessment of the imputation models is added.
