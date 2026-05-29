library(tidyverse)

# data----
brood <- read.csv("data/Situk_sockeye.csv", header = TRUE) %>% 
  dplyr::select(S = spawn, R = recruit50) %>% 
  mutate(n()) -> sr

nyrs <- nrow(sr) #calculates the number of years of data

# determine the dataset to run (only run one dat)
# AR1 Ricker
#dat <- list(nyrs = nyrs, S = sr$S, R = sr$R1, 
            #rk =1, ar1=1) 

# Basic Ricker
dat <- list(nyrs = nyrs, S = sr$S, R = sr$R, rk = 1, ar1 = 0, d = 4) 
# if running basic Ricker model, make sure that correct lnalpha.c chosen in the model.txt file
# ar1 = 0 and kf =0 in standard analysis   
# ar1 = 1 and kf =0 AR1 error moddel 
# ar1 = 0 and kf =1 time-varying alpha