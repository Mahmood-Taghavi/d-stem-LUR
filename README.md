# d-stem-LUR
This is the companion repository for the paper titled "Concurrent Spatiotemporal Daily Land Use Regression Modeling and Missing Data Imputation of Fine Particulate Matter Using Distributed Space-Time Expectation Maximization". 
Currently, we have added the dataset and cleaned code to fit and evaluate the final D-STEM kriging and D-STEM LUR models to the repository. We are going to add estimated pollution maps in upcoming updates. 
## data
This folder contains daily PM2.5 concentrations in 2015 and the corresponding pool of 210 potentially predictive variables (PPVs) in 30 monitoring stations in Tehran
## d-stem
This folder contains the original d-stem software source code which is required for model building
## code
This folder contains scripts and function we are coded to fit D-STEM spatiotemporal kriging and LUR models and to assess them using various metrics (19 metrics) and from different aspects (spatiotemporal, spatial, and temporal) in h-block cross-validation
