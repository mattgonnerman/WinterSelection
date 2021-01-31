### Data Management ###
## Takes a while, outputs saved so no need to run multiple times.
## First need to have a behavioral state to associate with each point
source(file = "Winter Selection - 1 HMM Forage Loaf Roost.R")

## Next need all the appropriate covariates of interest
source(file = "Winter Selection - 2 Data Management.R")

### Run Analysis ###
##Perform the INLA Analysis for SSF
source(file = "Winter Selection - 3 INLA Analysis.R")

### Results ###
## Create Graphs
source(file = "Winter Selection - 4 Graphs.R")