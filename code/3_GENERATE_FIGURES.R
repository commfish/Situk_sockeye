# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)
lower_bounds <- 33200 #lower bound of recommended escapement goal range (update)
upper_bounds <- 59300 #upper bound of recommended escapement goal range (update)

stat_quants <- read.csv(file= paste0(out.path,"/statsquants.csv"))
quant_lambert <- read.csv(file= paste0(out.path,"/quantiles_lambert.csv")) 

lnalpha <- stat_quants[stat_quants$variable == "lnalpha", "X50."]
lnalpha.c <- stat_quants[stat_quants$variable == "lnalpha.c", "X50."]
beta <- stat_quants[stat_quants$variable == "beta", "X50."]
beta1 <- beta*10^-4

SMSY <- quant_lambert[quant_lambert$variable == "Smsy_lambert", "X50"]
UMSY <- quant_lambert[quant_lambert$variable == "Umsy_lambert", "X50"]
SMAX <- quant_lambert[quant_lambert$variable == "Smax", "X50"] 
SEQ  <- quant_lambert[quant_lambert$variable == "Seq", "X50"]

# load----
library(plyr)
library(tidyverse)
library(reshape2)
library(grid)
library(gsl)
library(scales)
library(backports)
devtools::install_github("commfish/fngr")
library(fngr)
# font_import() # only need run once
# extrafont::font_import()
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_report(base_size = 14))
source('code/functions.r')

if(!dir.exists(file.path("output", "base_case", "processed"))){dir.create(file.path("output", "base_case", "processed"))}

# data----
read.csv("data/Situk_sockeye.csv") %>%
  mutate(yield = (recruit50 - spawn)) -> spawnrecruitdat
read.csv(file = paste0(out.path,"/coda.csv")) -> coda
read.csv("output/base_case/processed/recruit_data.csv") -> recruit

# analysis----
# function for probability profiles and figures
#profile(i = 10, z = 50, xa.start = 0, xa.end = 8000, lnalpha.c, beta1, coda) # can change i,z, xa.start, xa.end 
profile(i = 5, z = 20, xa.start = 0, xa.end = 20000, lnalpha.c, lnalpha, beta1, coda) # can change i,z, xa.start, xa.end 
QM <- read.csv("output/base_case/processed/QM.csv")
CI_median <- read.csv("output/base_case/processed/CI_median.csv")
CI_mean <- read.csv("output/base_case/processed/CI_mean.csv")

# horsetail (spawner recruit) plots
ggplot(data = CI_median, aes(x = Escapement, y = Median)) +
  geom_line(size=0.75, lty=1) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.1) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.2) +
  geom_line(data = CI_mean, aes(x = Escapement, y = Median), size=0.75, lty=2) + 
  xlab("Spawners (S)") +
  ylab("Recruits (R)") +
  #geom_vline(xintercept = SMSY, color ="gray70", lty=2) +
  scale_y_continuous(labels = comma,breaks = seq(0, 300000, 50000), limits = c(0, 300000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 250000, 50000), limits = c(0, 250000)) +
  geom_line(data = CI_median, aes(x = Escapement, y =Escapement),linetype="solid", size=0.75, color ="grey60") +
  geom_point(data = spawnrecruitdat, aes(x=spawn, y=recruit50),pch=1, size=2) + 
  theme(text=element_text(size=14)) +
  geom_text(size=3, data=spawnrecruitdat, aes(x=spawn, y=recruit50, label = year,family="Times",
                                             hjust = -0.1, vjust= -0.4))
ggsave("output/base_case/processed/horsetail.png", dpi = 500, height = 6, width = 8, units = "in")

# expected yield plot
ggplot(QM, aes(Escapement, Median))+geom_line(size=1)+
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.15)+
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.15)+ xlab("Escapement")+
  ylab("Expected Yield")+scale_y_continuous(labels = comma)+
  scale_x_continuous(labels = comma,breaks = seq(0, 250000, 50000), limits = c(0,250000))+
  scale_y_continuous(labels = comma,breaks = seq(-200000, 200000, 50000), limits = c(-200000,200000))+
  geom_hline(yintercept = 0, lty = 2, col = "grey10") + geom_point(data=spawnrecruitdat, aes(x=spawn, y=yield),pch=16, size = 2) 
ggsave("output/base_case/processed/expected_yield.png", dpi = 500, height = 6, width = 8, units = "in")

# yield profiles 
read.csv("output/base_case/processed/optimal_yield_data.csv") %>%
  filter(max_pct == 0.9) -> optimal_yield_data

optimal_yield_data %>%
  filter(sra == "Yield Profile") -> fig_data1

optimal_yield_data %>%
  filter(sra == "Overfishing Profile") -> fig_data2

optimal_yield_data %>%
  filter(sra == "Recruitment Profile") -> fig_data3

ggplot(fig_data1, aes(x = Escapement, y = Probability)) + ggtitle("(c) Yield Profile") + 
  annotate("rect", xmin = lower_bounds, xmax = upper_bounds, ymin = 0, ymax = 1,
           inherit.aes = FALSE, fill = "grey80", alpha = 0.9) +
  geom_line() +    theme(plot.title = element_text(size = 12, face = "bold"),
                         strip.text.y = element_text(size=0),legend.position= "none") +
  scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1))+
  scale_linetype_discrete(name = "Percent of Max.") + xlab('Escapement (S)')+
  facet_grid(sra ~ .)  -> plot1

ggplot(fig_data2, aes(x = Escapement, y = Probability)) + 
  annotate("rect", xmin = lower_bounds, xmax = upper_bounds, ymin = 0, ymax = 1,
           inherit.aes = FALSE, fill = "grey80", alpha = 0.9) + ggtitle("(a) Overfishing Profile") + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size=0),legend.position=c(0.92,0.86), legend.title = element_blank()) +
  geom_line() + xlab('Escapement (S)') +
  scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
  scale_linetype_discrete(name = "Percent of Max.") + 
  facet_grid(sra ~ .) -> plot2

ggplot(fig_data3, aes(x = Escapement, y = Probability)) + 
  annotate("rect", xmin = lower_bounds, xmax = upper_bounds, ymin = 0, ymax = 1,
           inherit.aes = FALSE, fill = "grey80", alpha = 0.9) + ggtitle("(b) Recruitment Profile") + 
  geom_line() + xlab('Escapement (S)') +   theme(plot.title = element_text(size = 12, face = "bold"),
                                                 strip.text.y = element_text(size=0),legend.position= "none") +
  scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
  scale_linetype_discrete(name = "Percent of Max.") +
  facet_grid(sra ~ .)  -> plot3
cowplot::plot_grid(plot2,plot3,plot1, align = "v", nrow = 3, ncol=1) 
ggsave("output/base_case/processed/yield_profiles.png", dpi = 500, height = 7, width = 6, units = "in")

# residual plot 
tickryr <- data.frame(Year = 1976:2020)
xaxis <- fngr::tickr(tickryr, Year, 4)
recruit <- read.csv("output/base_case/processed/recruit_data.csv") %>%
  dplyr::select(RD2.1:RD2.43)

names(recruit) <- as.character(1976:2018)

resid_summary <- recruit %>%
  pivot_longer(cols = everything(), names_to = "year", values_to = "resid") %>%
  mutate(year = as.integer(year)) %>%
  group_by(year) %>%
  dplyr::summarise(
    average = mean(resid, na.rm = TRUE),
    cil     = quantile(resid, 0.025, na.rm = TRUE),
    ciu     = quantile(resid, 0.975, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(year = as.numeric(year)) %>%
ggplot(., aes(year, average))+geom_line(linewidth=1) + geom_point(size=3)+
  geom_ribbon(aes(ymin = cil, ymax = ciu), alpha=.15)+ xlab("Year")+ylab("Residuals")+
  scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1976, 2020))+
  annotate("text",x = 1976, y=1.96, label="b)", family="Times" ,size=6) +theme (axis.text.x=element_text (size=10)) +
  scale_y_continuous(labels = comma,breaks = seq(-2, 2, 0.5), limits = c(-2,2))+
  geom_hline(yintercept = 0, color = "grey70", lty = 2) -> plot1

# predicted plot (recruits)
recruit <- read.csv("output/base_case/processed/recruit_data.csv") %>%
  dplyr::select(RD2.1:RD2.43)

names(recruit) <- as.character(1976:2018)

recruit %>%
  pivot_longer(cols = everything(), names_to = "year", values_to = "resid") %>%
  mutate(year = as.integer(year)) %>%
  group_by(year) %>%
  dplyr::summarise(
    average = mean(resid, na.rm = TRUE),
    cil     = quantile(resid, 0.025, na.rm = TRUE),
    ciu     = quantile(resid, 0.975, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(year = as.numeric(year)) %>%
  merge(., spawnrecruitdat, by="year") %>%
  mutate(
    logR  = log(recruit50),
    R.m   = logR - average,
    upper = logR - ciu,
    lower = logR - cil) %>%
ggplot(aes(year, R.m)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15) +
  geom_point(aes(x = year, y = logR), pch = 1, size = 3) +
  xlab("Year") + ylab("ln(Recruits)") +
  scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1976, 2020)) +
  scale_y_continuous(labels = comma, breaks = seq(10, 15, 0.5), limits = c(10, 15)) +
  annotate("text", x = 1976, y = 14.8, label = "a)", family = "Times", size = 6) +
  theme(axis.text.x = element_text(size = 10)) -> plot2
cowplot::plot_grid(plot2, plot1,  align = "v", nrow = 2, ncol=1)
ggsave("output/base_case/processed/resids_recruit.png", dpi = 500, height = 7, width = 8, units = "in")

# create optimal yield plot
read.csv("output/base_case/processed/yield.csv") -> df
get_two_nearest_bounds(df, "oy_0.9", 0.70)
get_two_nearest_bounds(df, "oy_0.9", 0.75)
get_two_nearest_bounds(df, "oy_0.9", 0.80)
get_two_nearest_bounds(df, "oy_0.9", 0.85)
get_two_nearest_bounds(df, "oy_0.9", 0.90)                     


