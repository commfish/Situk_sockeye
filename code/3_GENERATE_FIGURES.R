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
spawnrecruitdat <- read.csv("data/Situk_sockeye.csv") 
coda <- read.csv(file = paste0(out.path,"/coda.csv"))

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
  scale_y_continuous(labels = comma,breaks = seq(0, 200000, 50000), limits = c(0, 200000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
  ylab("Recruits (R)")+xlab("Spawners (S)")+
  geom_line(data=dataset, aes(x=Escapement, y=Escapement, group=1),linetype="solid", size=1)+ # replacement line
  geom_line(data=subset(dataset,variable==51),colour = "black", lty=2, size=1)+
  geom_point(data=subset(dataset,variable==52),colour = "black", pch=16, size=1)+
  geom_text(size=3, data=dataset, aes(x=Escapement, y=recruitment, group=52, label=year,family="Times", 
                                      hjust = -0.1, vjust= -0.4)) 
ggsave("output/base_case/processed/horsetail.png", dpi = 500, height = 6, width = 8, units = "in")

# horsetail (spawner recruit) plots
ggplot(data = CI, aes(x = Escapement, y = Median)) +
  geom_line(size=0.75, lty=2) +
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.1) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.2) +
  xlab("Spawners (S)") +
  ylab("Recruits (R)") +
  geom_vline(xintercept = SMSY, color ="gray70", lty=2) +
  scale_y_continuous(labels = comma,breaks = seq(0, 200000, 50000), limits = c(0, 200000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 400000, 50000), limits = c(0, 400000)) +
  geom_line(aes(x = Escapement, y =Escapement),linetype="solid", size=0.75, color ="black") +
  geom_point(data = spawnrecruitdat, aes(x=spawn, y=recruit),pch=1, size=2) + 
  theme(text=element_text(size=14)) +
  geom_text(size=3, data=spawnrecruitdat, aes(x=spawn, y=recruit, label=year,family="Times",
                                             hjust = -0.1, vjust= -0.4)) -> plot1
cowplot::plot_grid(plot1,  align = "v", nrow = 1, ncol=1)
ggsave("output/base_case/processed/horsetail_alt.png", dpi = 500, height = 6, width = 8, units = "in")