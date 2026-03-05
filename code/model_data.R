library(tidyverse)

# data----
brood <- read.csv("data/Situk_sockeye.csv", header = TRUE)
brood

# cleanup data
brood %>% 
  dplyr::select(S = spawn, R = recruit) %>% 
  mutate(n()) -> sr

nyrs <- nrow(sr) #calculates the number of years of data.
dat <- list(nyrs = nyrs, S = sr$S, R = sr$R)
