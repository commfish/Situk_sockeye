# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)

LowerB <- 38000 #lower bound of recommended escapement goal range
UpperB <- 86000 #upper bound of recommended escapement goal range
SMSY <- 52774.07 #Lambert W version of SMSY from file
lnalpha.c <- 2.143810675
lnalpha <- 1.685389252
beta <- 0.000014123254157  #let's try and grab these last 4 parameter outputs from the "stats.csv" output

#load----
library(plyr)
library(tidyverse)
library(reshape2)
library(grid)
library(gsl)
library(scales)
library(backports)
devtools::install_github("commfish/fngr")
library(fngr)
#library(FField)
#font_import() # only need run once
extrafont::font_import()
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_report(base_size = 14))
source('code/functions.r')

if(!dir.exists(file.path("output", "base_case", "processed"))){dir.create(file.path("output", "base_case", "processed"))}

# data----
spawnrecruitdat <- read.csv("data/Situk_sockeye.csv") #Load Data File (make sure this file is updated)
coda <- read.csv("output/base_case/coda.csv") 

# data clean----
# profile spawnrecruitdat
coda %>% 
  mutate(S.eq.c = lnalpha.c/beta, 
         S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
         MSY.c = R.msy.c-S.msy.c, 
         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda

# analysis----
# create function for probability profiles and figures
profile(i=10, z=50, xa.start=0, xa.end=8000,lnalpha.c, beta) #can change i,z, xa.start, xa.end
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
          year = "")%>%
  rbind(., dataset) -> dataset

spawnrecruitdat %>%
  mutate(Escapement = spawn,
         lnalpha.c = NA,
         beta = NA,
         recruitment =recruit,
         variable = 52) %>%
  dplyr::select(c(Escapement, lnalpha.c, beta, recruitment, variable, year)) %>%
  rbind(., dataset) %>%
  arrange(variable, Escapement) -> dataset

ggplot(data=dataset, aes(x=Escapement, y=recruitment, group=variable))+
  geom_line(data=subset(dataset,dataset$variable<52),linetype="solid", size=0.5, color="grey80")+
  scale_y_continuous(labels = comma,breaks = seq(0, 1250000, 250000), limits = c(0, 1250000))+
  scale_x_continuous(labels= comma,breaks = seq(0, 400000, 100000), limits = c(0, 400000))+
  ylab("Recruits (R)")+xlab("Spawners (S)")+
  geom_line(data=dataset, aes(x=Escapement, y=Escapement, group=1),linetype="solid", size=1)+#replacement line
  geom_line(data=subset(dataset,variable==51),colour = "black", lty=2, size=2)+
  geom_point(data=subset(dataset,variable==52),colour = "black", pch=16, size=1)+
  geom_text(size=3, data=dataset, aes(x=Escapement, y=recruitment, group=52, label=year,family="Times", 
                                      hjust = -0.1, vjust= -0.4)) 
ggsave("output/base_case/processed/horsetail.png", dpi = 500, height = 6, width = 8, units = "in")

# horsetail (spawner recruit) plots

ggplot(data=CI, aes(x=Escapement, y=Median)) +
  geom_line(size=0.75, lty=2) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.05) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.2) +
  xlab("Spawners (S)") +
  ylab("Recruits (R)") +
  geom_vline(xintercept = SMSY, color ="gray70", lty=2)+
  scale_y_continuous(labels = comma,breaks = seq(0, 500000, 50000), limits = c(0, 500000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 350000, 50000), limits = c(0, 350000)) +
  geom_line(aes(x=Escapement, y=Escapement),linetype="solid", size=0.75, color ="black") +
  geom_point(data=spawnrecruitdat, aes(x=spawn, y=recruit),pch=1, size=2) + # need this data!
  theme(text=element_text(size=17)) +
  geom_text(size=3, data=spawnrecruitdat, aes(x=spawn, y=recruit, label=year,family="Times",
                                             hjust = -0.1, vjust= -0.4)) -> plot1
cowplot::plot_grid(plot1,  align = "v", nrow = 1, ncol=1)
ggsave("output/base_case/processed/horsetail_alt.png", dpi = 500, height = 6, width = 8, units = "in")