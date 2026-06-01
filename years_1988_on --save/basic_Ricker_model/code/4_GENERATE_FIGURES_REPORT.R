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
#pak::pak("adamreimer/EGprocess")
library(EGprocess)
devtools::install_github("commfish/fngr")
library(fngr)
# font_import() # only need run once
# extrafont::font_import()
# data
read.csv("data/Situk_sockeye_run.csv") -> Situk_sockeye_run # this dataset includes 1976-1987
read.csv("data/Situk_sockeye_run_historic.csv") -> Situk_sockeye_run_historic 
# read.csv("data/Situk_sockeye.csv") %>%
#   mutate(S = spawn,
#          R = recruit50,
#          yr = year,         
#          Y = (recruit50 - spawn)) %>%
#   dplyr::select(yr, S, R, Y) -> Situk_sockeye_76

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

read.csv("years_1988_on/basic_Ricker_model/output/post_data.csv") -> post_Situk_byr88_18 
# read.csv("years_1976_on/basic_Ricker_model/output/post_data.csv") -> post_Situk_byr76_18

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
                                          escapement goal are indicated with black fill.", width = 85))
out.file <- paste0(out.path, "/output/processed/Esc.png")
ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")
# spawner recruit plot (see table 12 in Clark et al. 2002); sigma is square root of the residual mean square error
# 1976-2018
# post_Situk_byr76_18 <- list("Brood years: 1976-2018" = post_Situk_byr76_18)
# 
# post_onlyparameters <-
#   c(list('Brood years: 1976-1997' = c(lnalpha = 1.396244692, beta = 0.11, sigma = 0.36)),
#     post_Situk_byr76_18)
# 
# goal_onlyparameters_new <-goal_Situk %>%
#   bind_rows(data.frame(yr = 2026,
#                        lb = 27500,
#                        ub = 70000))
# 
# plot_SR <- plot_SR(post_onlyparameters,
#         Situk_sockeye_historic,
#         goal_dat = goal_onlyparameters_new,
#         "Situk River Sockeye Salmon",
#         new_finding = TRUE,
#         multiplier = 1e-4)
# plot_SR + geom_point(data = Situk_sockeye_76, pch =16, size =2) +
#   labs(subtitle = paste0("Brood Years: 1976-2018"))
# out.file <- paste0(out.path, "/output/processed/SR_76_18.png")
# ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")
# spawner recruit plot (see table 12 in Clark et al. 2002); sigma is square root of the residual mean square error
# 1988-2018
post_Situk_byr88_18 <- list("Brood years: 1988-2018" = post_Situk_byr88_18)

post_onlyparameters <-
  c(list('Brood years: 1976-1997' = c(lnalpha = 1.396244692, beta = 0.11, sigma = 0.36)),
    post_Situk_byr88_18)

goal_onlyparameters_new <-goal_Situk %>%
  bind_rows(data.frame(yr = 2026,
                       lb = 27500,
                       ub = 70000))

plot_SR <- plot_SR(post_onlyparameters,
                   Situk_sockeye_historic,
                   goal_dat = goal_onlyparameters_new,
                   "Situk River Sockeye Salmon",
                   new_finding = TRUE,
                   multiplier = 1e-4)
plot_SR + geom_point(data = Situk_sockeye_88, pch =16, size =2) +
  labs(subtitle = paste0("Brood Years: 1988-2018"))

out.file <- paste0(out.path, "/output/processed/SR_88_18.png")
ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")
       
# expected yield plot
posterior_missing <- 
  c(list('Brood years: 1976-1997' = NULL),
    post_Situk_byr88_18)

goal_update <-
  goal_Situk %>%
  bind_rows(
    data.frame(yr = 2026,
               lb = 27500,
               ub = 70000))


plot_EY <-plot_ey(posterior_missing,
        Situk_sockeye_historic,
        goal_data = goal_update,
        "Situk River Sockeye Salmon",
        new_finding =TRUE,
        multiplier = 1e-4)
plot_EY + geom_point(data = Situk_sockeye_88, pch =16, size = 2) +
  labs(subtitle = paste0("Brood Years: 1988-2018"))+ labs(caption = stringr::str_wrap(
    "Note: Hollow circles indicate the data when the escapement goal last changed, while filled circles 
    and solid lines indicate the data collected since and the estimate of median sustained yield from all available data. Vertical lines show the escapement
    that maximizes sustained yield. The new escapement goal finding is shaded brown.", width = 85))
out.file <- paste0(out.path, "/output/processed/EY.png")
ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")

# optimal yield plot
profile_missing <- get_profile(post_Situk_byr88_18, multiplier = 1e-4)

plot_profile(profile_missing,
             goal_update,
             new_finding =TRUE,
             "Situk River Sockeye Salmon")
out.file <- paste0(out.path, "/output/processed/profile.png")
ggsave(out.file, dpi = 500, height = 6, width = 8, units = "in")

# spawner recruit table
table_SR(post_Situk_byr88_18[2], 
         title = "Situk River Sockeye Salmon", 
         multiplier = 1e-4)