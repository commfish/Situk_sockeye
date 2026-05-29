# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)
# load----
library(plyr)
library(tidyverse)
library(reshape2)
library(grid)
library(gsl)
library(scales)
library(backports)
library(ggplot2)
#install.packages("pak")
pak::pak("adamreimer/EGprocess")
library(EGprocess)
devtools::install_github("commfish/fngr")
library(fngr)
# font_import() # only need run once
# extrafont::font_import()
# data
read.csv("data/Situk_sockeye_run.csv") -> Situk_sockeye_run
read.csv("data/Situk_sockeye_run_historic.csv") -> Situk_sockeye_run_historic 
read.csv("data/Situk_sockeye.csv") %>%
  mutate(S = spawn,
         R = recruit50,
         yr = year,         
         Y = (recruit50 - spawn)) %>%
  dplyr::select(yr, S, R, Y) -> Situk_sockeye_76

read.csv("data/Situk_sockeye.csv") %>%
  mutate(year = as.numeric(year)) %>%
  filter(year > 1987) %>%
  mutate(S = spawn,
         R = recruit50,
         yr = year,         
         Y = (recruit50 - spawn)) %>%
  dplyr::select(yr, S, R, Y)  -> Situk_sockeye_88

read.csv("data/Situk_sockeye_historic.csv") %>%
  mutate(S = spawn,
         R = recruit50,
         yr = year,
         Y = (recruit50 - spawn)) %>%
  dplyr::select(yr, S, R, Y) -> Situk_sockeye_historic 

read.csv("data/goal_Situk.csv") -> goal_Situk 


# create brood tables
brood_Situk <- get_brood(run_data = Situk_sockeye_run)
brood_Situk_historic <- get_brood(run_data = Situk_sockeye_run_historic)

# plot escapement
plot_S <- plot_escapement(brood_data = brood_Situk, goal_data = goal_Situk, 
                      title = "Situk River Sockeye Salmon")
plot_S + labs(caption = stringr::str_wrap("Note: Escapement goal lower and 
                                          upper bounds are shown as solid and dashed
                                          lines, respectively. Escapements below
                                          the lower bound of the contemporaneous
                                          escpements goal are indicated with black fill. 
                                          Goals prior to 1995...", width = 85))

# spawner recruit plot (see table 12 in Clark et al. 2002); sigma is square root of the residual mean square error
read.csv(file = paste0(out.path,"/output/post_data.csv")) -> post_Situk_byr76_18

post_Situk_byr76_18 <- list("Brood years: 1976-2018" = post_Situk_byr76_18)

post_onlyparameters <-
  c(list('Brood years: 1976-1997' = c(lnalpha = 1.396244692, beta = 0.11, sigma = 0.36)),
    post_Situk_byr76_18)

goal_onlyparameters_new <-goal_Situk %>%
  bind_rows(data.frame(yr = 2026,
               lb = 25000,
               ub = 100000))

plot_SR <- plot_SR(post_onlyparameters,
        Situk_sockeye_historic,
        goal_dat = goal_onlyparameters_new,
        "Situk River Sockeye Salmon",
        new_finding = TRUE,
        multiplier = 1e-4)
plot_SR + geom_point(data = Situk_sockeye, pch =16, size =2) +
  labs(subtitle = paste0("Brood Years: 1976-2018"))
       
# expected yield plot
posterior_missing <- 
  c(list('Brood years: 1976-1997' = NULL),
    post_Situk_byr76_18)

goal_new <- data.frame(yr = 2026, lb = 250000, ub = 100000)

plot_EY <-plot_ey(posterior_missing,
        Situk_sockeye_historic,
        goal_dat = goal_Situk,
        "Situk River Sockeye Salmon",
        multiplier = 1e-4)
plot_EY + geom_point(data = Situk_sockeye, pch =16, size = 2) +
  labs(subtitle = paste0("Brood Years: 1976-2018"))

# optimal yield plot
profile_missing <- get_profile(post_Situk_byr76_18, multiplier = 1e-4)

plot_profile(profile_missing,
             goal_Situk,
             "Situk River Sockeye Salmon")

# spawner recruit table
table_SR(post_Situk_byr76_18[2], 
         title = "Situk River Sockeye Salmon", 
         multiplier = 1e-4)