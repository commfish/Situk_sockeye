library(tidyverse)

# data----
brood <- read.csv("years_1976_on/basic_Ricker_model/data/Situk_sockeye.csv", header = TRUE)
brood

# cleanup data
brood %>% 
  dplyr::select(S = spawn, R1 = recruit25, R2 = recruit50, R3 = recruit75) %>% 
  mutate(n()) -> sr

nyrs <- nrow(sr) #calculates the number of years of data

# determine the dataset to run (only run one dat)
# AR1 Ricker
#dat <- list(nyrs = nyrs, S = sr$S, R = sr$R1, 
            #rk =1, ar1=1) 

# Basic Ricker
dat <- list(nyrs = nyrs, S = sr$S, R = sr$R2, rk = 1, ar1 = 0, d = 4) 
# if running basic Ricker model, make sure that correct lnalpha.c chosen in the model.txt file
