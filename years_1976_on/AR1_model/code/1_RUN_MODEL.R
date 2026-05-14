# Situk sockeye state space model
# authors: Sara E Miller 
# Last edited: March 2026
# must download program JAGS for this script to work

# load libraries----
library(coda)
library(tidyverse)
library(R2jags)
library(rjags)
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
library(runjags)
library(grid)
library(gridExtra)

# STEP 1: CHOOSE SETTINGS----

# if test runs then do sensitivity tests with explore, and final run with full
# "explore" version takes ~10min with the current settings.
package.use <- "rjags"  
jags.settings <- "test"  # "test" or "explore" or full" 
sensitivity.analysis <- 0 #0; 1 is yes and 0 is no

# load custom functions
source("years_1976_on/AR1_model/code/functions.R")

# create output folder for model results
out.path <- paste0("years_1976_on/AR1_model/output/")
if(!exists(out.path)){dir.create(out.path)}

# create output folder for model results
out.path <- paste0("years_1976_on/AR1_model/output/processed/")
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

if (jags.settings == "full") {
  n.adapt.use <- 20000      # number of adaptation iterations (tuning phase)
  n.iter.use  <- 200000     # total iterations per chain (after adaptation)
  n.burnin.use <- 20000     # burn-in iterations (discarded before sampling)
  thin.use    <- 100        # thinning interval (keep every 100th sample)
  by.use      <- 1000       # step size for progress bar updates
}


# STEP 2: READ IN DATA, MODEL, AND INITIAL VALUES----
# generates the object "dat"
source("years_1976_on/AR1_model/code/model_data.R")

# generates initial values
source("years_1976_on/AR1_model/code/model_inits.R")

# STEP 3: RUN THE MODEL AND PROCESS THE OUTPUT----
start.jags <- proc.time()
sw.randomseed <- 200
if(package.use == "rjags" & sensitivity.analysis == 0){
  parameters <- c("lnalpha", "phi", "beta", "sigma", "sigmaw", "Tau", "tauw", "alpha", "lnalpha.c", "e0", "e")
  
  jmod <- rjags::jags.model(
    file='years_1976_on/AR1_model/code/Situk_sockeye.txt', 
    data = dat, n.chains = 3, 
    inits = inits, 
    n.adapt = n.adapt.use) 
  
  stats::update(jmod, n.iter = n.iter.use, by = by.use, progress.bar = 'text', DIC=T, n.burnin = n.burnin.use) # this modifies the original object, function returns NULL
  post <- rjags::coda.samples(jmod, parameters, n.iter = n.iter.use, thin = thin.use, n.burnin = n.burnin.use)
  
  end.jags <- proc.time()   # store time for MCMC
  post.arr <- as.array(post) # convert to an accessible obj
  
# run the script that generates all the outputs 
  end.output  <- proc.time() 
  print("Jags took")
  print(end.jags - start.jags)
  print("Output Processing took")
  print(end.output - end.jags)
}



