# Situk sockeye state space model
# authors: Sara E Miller 
# Last edited: SMarch 2026
# must download program JAGS for this script to work

# load libraries----
library(coda)
library(tidyverse)
library(R2jags)
library(rjags)
library(runjags)
library(R2OpenBUGS)
library(gsl)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(gdata)
library(ggplot2)
library(scales)
library(cowplot)
library("devtools")
devtools::install_github("commfish/fngr")
library(fngr)
library(extrafont)

# STEP 1: CHOOSE SETTINGS----

# if test runs then do sensitivity tests with explore, and final run with full
# "explore" version takes ~10min with the current settings.
out.label <-  "base_case"  # label to be used for the output folder (and for scenario comparisons)
package.use <- "rjags"  #"rjags"  or "R2jags"
jags.settings <- "test"  # "test" or "explore" or full" 
sensitivity.analysis <- 0 #0; 1 is yes and 0 is no

source("code/model_source.R") 
print(jag.model.SR)
model_file_loc = paste("code/","Situk_sockeye.txt", sep="")
write.model(jag.model.SR, model_file_loc)

# load custom functions
source("code/functions.R")

# create output folder for model results
out.path <- paste0("output/", out.label)
if(!exists(out.path)){dir.create(out.path)}

# create output folder for model results
out.path <- paste0("output/", out.label)
if(!exists(out.path)){dir.create(out.path)}

# choices of model runs
if(jags.settings == "test"){
  n.adapt.use <- 100 ; n.iter.use <- 500;  n.burnin.use <- 100;   thin.use = 10
  by.use <- 10 # this is just for the progress bar
}

if(jags.settings == "explore"){
  n.adapt.use <- 10000 ; n.iter.use <- 10000; n.burnin.use <- 30000 ;   thin.use = 10
  by.use <- 100 # this is just for the progress bar
}

if(jags.settings == "full"){
  n.adapt.use <- 10000  ; n.iter.use <- 1000000    #1,000,000 per chain; 3 chains; thin by 1000
  n.burnin.use <- 100000  # consider increasing this?
  thin.use = 1000; by.use <- 1000 # this is just for the progress bar 
}

# STEP 2: READ IN DATA, MODEL, AND INITIAL VALUES----
# generates the object "dat"
source("code/model_data.R")

# generates initial values
source("code/model_inits.R")

# STEP 3: RUN THE MODEL AND PROCESS THE OUTPUT----
# rjags
start.jags <- proc.time()
sw.randomseed <- 2540
if(package.use == "rjags" & sensitivity.analysis == 0){
  parameters <- c("lnalpha", "phi", "beta", "sigma", "sigmaw", "Tau", "tauw", "alpha", "lnalpha.c")
  
  jmod <- rjags::jags.model(file='code/Situk_sockeye.txt', data = dat, n.chains = 3, inits = inits, n.adapt = n.adapt.use) 
  stats::update(jmod, n.iter = n.iter.use, by = by.use, progress.bar = 'text', DIC=T, n.burnin = n.burnin.use) # this modifies the original object, function returns NULL
  post <- rjags::coda.samples(jmod, parameters, n.iter = n.iter.use, thin = thin.use, n.burnin = n.burnin.use)
  
  end.jags <- proc.time()   # store time for MCMC
  post.arr <- as.array(post) # convert to an accessible obj
  
  # run the script that generates all the outputs 
  #source("state_space_model/2023_analysis/recommendations_added/code/2_GENERATE_OUTPUTS.R")
  end.output  <- proc.time() 
  print("Jags took")
  print(end.jags - start.jags)
  print("Output Processing took")
  print(end.output - end.jags)
}

