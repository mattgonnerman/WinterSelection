### Data Management ###
## Prep Movement Data and Covariates for momentuHMM
source(file = "Winter Selection - 1a HMM Data Prep.R")

## Run Candidate Models
source(file = "Winter Selection - 1b Run HMM Models.R")

## HMM Model Selection and Graphing Options
source(file = "Winter Selection - 1c HMM Results.R")

# source(file = "Winter Selection - 1d HMM Graphs.R")

# Run just the best HMM model (Once Identified)


## Next need all the appropriate covariates of interest
source(file = "Winter Selection - 2 SSF Data Management.R")

##Perform the INLA Analysis for SSF
source(file = "Winter Selection - 3 INLA Analysis.R")

## Organize Results in a more readable format
source(file = "Winter Selection - 4 Results Clean Up.R")

## Create Graphs
#Doesn't appear to run when source is used, but fine when run from script directly
source(file = "Winter Selection - 5 Graphs.R")
