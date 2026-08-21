# Situk Summary Figures
# authors: Justin Priest
# Last edited: August 2026


library(tidyverse)
library(EGprocess)
#remotes::install_github("commfish/adfggraph")
library(adfggraph)


rundata <- read_csv("data/Situk_sockeye_run.csv") %>%
  rename(Year = yr, escapement = S)

runharvdata <- read_csv("data/Situk_sockeye_harvestsources.csv") %>%
  select(Year, commharv_total, subharv_total, sportharv_total, allharv_total) %>%
  left_join(rundata %>% select(Year, escapement)) %>%
  mutate(totalrun = allharv_total + escapement)
runharvdata

runharvprop <- runharvdata %>% 
  mutate(commharvrate  = commharv_total  / totalrun,
         subharvrate   = subharv_total   / totalrun,
         sportharvrate = sportharv_total / totalrun,
         allharvrate   = allharv_total   / totalrun) %>% 
  select(Year, commharvrate, subharvrate, sportharvrate, allharvrate) 
runharvprop 


plt_harvestrate <-runharvprop %>% 
  pivot_longer(-Year, names_to = "harvestsource", values_to = "proportion") %>%
  filter(harvestsource != "allharvrate") %>% 
  # reorder so that they appear in plot with sport on top, comm on bottom
  #  then rename them 
  mutate(harvestsource = fct_relevel(harvestsource,"sportharvrate", 
                           "subharvrate","commharvrate"),
         # rename the values 
         harvestsource = fct_recode(harvestsource,
                             "Sport   "       = "sportharvrate",
                             "Subsistence   " = "subharvrate",
                             "Commercial   "  = "commharvrate")) %>%
  ggplot(aes(x=Year, y = proportion, fill = harvestsource)) +
  geom_col() + 
  scale_x_continuous(breaks = seq(1980, 2025, 10),
                     minor_breaks = seq(1985, 2025, 5),
                     guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0.01, 0.05))) + #reduce gap at bottom
  scale_fill_manual(values = c(
      "Commercial   "  = "#0072B2",
      "Subsistence   " = "#D55E00",
      "Sport   "       = "#009E73")) +
  # alt palette
  # scale_fill_manual(values = c(
  #   "Commercial   "  = "#6C88A4",
  #   "Subsistence   " = "#C07A62",
  #   "Sport   "       = "#8FA97A"
  # )) +

  labs(x = "", y = "Harvest Rate", fill = "Fishery") + 
  guides(fill = guide_legend(nrow = 1, title.position = "top",
                             byrow = TRUE)) +
  theme_adfg(font_family = "Arial", font_size = 12, box = TRUE,
             legend.position = "bottom",
             legend.title.align = 0.5) # centered
plt_harvestrate

ggsave(plt_harvestrate, file = "output/harvestratesources.png", dpi = 300, 
       width = 6, height = 4, units = "in")

