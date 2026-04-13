# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)
LowerB <- 38000 #lower bound of recommended escapement goal range (update)
UpperB <- 86000 #upper bound of recommended escapement goal range (update)

df <- read.csv(file= paste0(out.path,"/statsquants.csv"))
df1 <- read.csv(file= paste0(out.path,"/quantiles_lambert.csv")) 

lnalpha <- df[df$variable == "lnalpha", "X50."]
lnalpha.c <- df[df$variable == "lnalpha.c", "X50."]
beta <- df[df$variable == "beta", "X50."]

SMSY <- df1[df1$variable == "Smsy_lambert.c", "X50"]
UMSY <- df1[df1$variable == "Umsy_lambert", "X50"]
SMAX <- df1[df1$variable == "Smax", "X50"] 
SEQ <- df1[df1$variable == "Seq", "X50"]

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

# analysis----
# function for probability profiles and figures
profile(i=10, z=50, xa.start=0, xa.end=8000,lnalpha.c, beta) # can change i,z, xa.start, xa.end 
QM <- read.csv("output/base_case/processed/QM.csv")
CI <- read.csv("output/base_case/processed/CI.csv")
num <- nrow(QM)
QM %>%
  dplyr::select(c(Escapement)) -> x
coda %>%
  dplyr::select(c(lnalpha.c, beta)) %>%
  filter(row_number()==1:50) %>%
  slice(rep(1:n(), each = num)) %>%
  cbind(., x) %>% #lnalpha.c, beta, and S 
  mutate(recruitment = Escapement*exp(lnalpha.c-beta*Escapement),
         variable = rep(1:50,each=num),
         year = "") -> dataset

QM %>%
  dplyr::select(c(Escapement)) %>%
  mutate (lnalpha.c = lnalpha.c,
          beta = beta,
          recruitment = Escapement*exp(lnalpha.c-beta*Escapement),
          variable = 51,
          year = "") %>%
  rbind(., dataset) -> dataset

spawnrecruitdat %>%
  mutate(Escapement = spawn,
         lnalpha.c = NA,
         beta = NA,
         recruitment =recruit50,
         variable = 52) %>%
  dplyr::select(c(Escapement, lnalpha.c, beta, recruitment, variable, year)) %>%
  rbind(., dataset) %>%
  arrange(variable, Escapement) -> dataset

# horsetail (spawner recruit) plots
ggplot(data = CI, aes(x = Escapement, y = Median)) +
  geom_line(size=0.75, lty=2) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.1) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.2) +
  xlab("Spawners (S)") +
  ylab("Recruits (R)") +
  geom_vline(xintercept = SMSY, color ="gray70", lty=2) +
  scale_y_continuous(labels = comma,breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
  geom_line(aes(x = Escapement, y =Escapement),linetype="solid", size=0.75, color ="grey60") +
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
  scale_x_continuous(labels = comma,breaks = seq(0, 200000, 50000), limits = c(0,200000))+
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
  annotate("rect", xmin = 40000, xmax = 75000, ymin = 0, ymax = 1,
           inherit.aes = FALSE, fill = "grey80", alpha = 0.9) +
  geom_line() +    theme(plot.title = element_text(size = 12, face = "bold"),
                         strip.text.y = element_text(size=0),legend.position= "none") +
  scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
  scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1))+
  scale_linetype_discrete(name = "Percent of Max.") + xlab('Escapement (S)')+
  facet_grid(sra ~ .)  -> plot1

ggplot(fig_data2, aes(x = Escapement, y = Probability)) + 
  annotate("rect", xmin = 40000, xmax = 75000, ymin = 0, ymax = 1,
           inherit.aes = FALSE, fill = "grey80", alpha = 0.9) + ggtitle("(a) Overfishing Profile") + 
  theme(plot.title = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(size=0),legend.position=c(0.92,0.86), legend.title = element_blank()) +
  geom_line() + xlab('Escapement (S)') +
  scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
  scale_linetype_discrete(name = "Percent of Max.") + 
  facet_grid(sra ~ .) -> plot2

ggplot(fig_data3, aes(x = Escapement, y = Probability)) + 
  annotate("rect", xmin = 40000, xmax = 75000, ymin = 0, ymax = 1,
           inherit.aes = FALSE, fill = "grey80", alpha = 0.9) + ggtitle("(b) Recruitment Profile") + 
  geom_line() + xlab('Escapement (S)') +   theme(plot.title = element_text(size = 12, face = "bold"),
                                                 strip.text.y = element_text(size=0),legend.position= "none") +
  scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
  scale_linetype_discrete(name = "Percent of Max.") +
  facet_grid(sra ~ .)  -> plot3
cowplot::plot_grid(plot2,plot3,plot1, align = "v", nrow = 3, ncol=1) 
ggsave("output/base_case/processed/yield_profiles.png", dpi = 500, height = 7, width = 6, units = "in")

# residual plot *****FIX DOWN BELOW
year <- data.frame(spawnrecruitdat$year)
data_Main %>%
  dplyr::select(RD2.1:RD2.32) -> data_Main
v1<- c(1983:2014)
names(data_Main)[1:32] <- v1

data_Main %>%
  gather(key = year, value = resid) %>%
  as.data.frame()%>%
  group_by(year) %>%
  summarise(average = mean(resid, na.rm=T),
            cil = quantile(resid, 0.025, na.rm=T),
            ciu = quantile(resid, 0.975, na.rm=T)) %>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(year = as.numeric(year)) %>%
  ggplot(., aes(year, average))+geom_line(linewidth=1) + geom_point(size=3)+
  geom_ribbon(aes(ymin = cil, ymax = ciu), alpha=.15)+ xlab("Year")+ylab("Residuals")+
  scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1982, 2016))+
  annotate("text",x = 1982, y=1.96, label="a)", family="Times" ,size=6) +theme (axis.text.x=element_text (size=10)) +
  scale_y_continuous(labels = comma,breaks = seq(-2, 2, 0.5), limits = c(-2,2))+
  geom_hline(yintercept = 0) -> plot1

# predicted plot (recruits)
read.csv(file = paste0("FIGURES_REPORT/output/processed/recruit_data_Main.csv")) -> data_Main
year <- data.frame(SR_Main_raw$Year)
data_Main %>%
  dplyr::select(RD.1:RD.32) -> data_Main
v2<- c(1983:2014)
names(data_Main)[1:32] <- v2

data_Main %>%
  gather(key = year, value = resid) %>%
  as.data.frame()%>%
  group_by(year) %>%
  summarise(average = mean(resid, na.rm=T),
            cil = quantile(resid, 0.025, na.rm=T),
            ciu = quantile(resid, 0.975, na.rm=T))%>%
  ungroup() %>%
  as.data.frame() %>%
  mutate(Year = as.numeric(year)) %>%
  merge(., dat, by="Year") %>%
  mutate(R.m = logR - average,
         upper = logR-cil,
         lower = logR-ciu) %>%
  ggplot(., aes(Year, R.m))+geom_line(linewidth=1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha=.15)+ xlab("Year")+ylab("ln(Recruits)")+
  scale_x_continuous(breaks = xaxis$breaks, labels = xaxis$labels, limits = c(1982, 2016))+
  annotate("text",x = 1982, y=11.97, label="b)", family="Times" ,size=6) + theme (axis.text.x=element_text (size=10)) +
  scale_y_continuous(labels = comma,breaks = seq(9, 12, 0.5), limits = c(9,12))+
  geom_point(aes(x=Year, y=logR),pch=1, size=3) -> plot2
cowplot::plot_grid(plot1, plot2,  align = "v", nrow = 2, ncol=1)
ggsave("FIGURES_REPORT/output/processed/resids_Main.png", dpi = 500, height = 7, width = 8, units = "in")

