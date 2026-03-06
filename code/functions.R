# https://rdrr.io/github/seananderson/ggsidekick/src/R/theme_sleek.R
theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}

# profile function----
# this function created by Ben William (Ben.Williams@alaska.gov) and adapted by Sara Miller (Sara.Miller@alaska.gov)
profile <-function(i,z,xa.start, xa.end,lnalpha.c, beta){
  xa = seq(xa.start, xa.end, by=i) 
  x =(xa+i)*z
  # empty dataframes
  dat <- data.frame(S0=rep(1, length(coda[,1])))
  dat1 <- data.frame(S0=rep(0, length(coda[,1])))
  dat2 <- data.frame(S0=rep(0, length(coda[,1])))
  dat3 <- data.frame(S0=rep(1, length(coda[,1])))
  dat4 <- data.frame(S0=rep(0, length(coda[,1])))
  dat5 <- data.frame(S0=rep(0, length(coda[,1])))
  dat6 <- data.frame(S0=rep(1, length(coda[,1])))
  dat7 <- data.frame(S0=rep(0, length(coda[,1])))
  dat8 <- data.frame(S0=rep(0, length(coda[,1])))
  dat9 <- data.frame(S0=rep(0, length(coda[,1])))
  dat10 <- data.frame(S0=rep(0, length(coda[,1])))
  for (i in 1:length(xa)){
    dat[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.7*coda$MSY.c), 0, ifelse(dat[,i]==0, 0,1))
    dat1[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.7*coda$MSY.c), 1,0)
    dat2[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.7*coda$Rmax), 1,0)
    dat3[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.8*coda$MSY.c), 0, ifelse(dat3[,i]==0, 0,1))
    dat4[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.8*coda$MSY.c), 1,0)
    dat5[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.8*coda$Rmax), 1,0)
    dat6[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.9*coda$MSY.c), 0, ifelse(dat6[,i]==0, 0,1))
    dat7[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i])-x[i])>(0.9*coda$MSY.c), 1,0)
    dat8[,i+1] = ifelse((x[i] * exp(coda$lnalpha.c-coda$beta*x[i]))>(0.9*coda$Rmax), 1,0)
    dat9[,i+1] = x[i]*exp(coda$lnalpha.c-coda$beta*x[i])-x[i] #expected yield
    dat10[,i+1] = x[i]*exp(coda$lnalpha.c-coda$beta*x[i]) # CI around S
  }
  # Overfishing estimate ----
  f.over <- function(x){
    x %>% 
      dplyr::filter(complete.cases(.)) %>% 
      dplyr::summarise_all(funs(mean)) %>% 
      tidyr::gather() %>% 
      dplyr::select(value)
  }
  
  of_0.7 <- f.over(dat)
  of_0.8 <- f.over(dat3)
  of_0.9 <- f.over(dat6)
  
  # Optimal yield estimate ----
  oy_0.7 <- f.over(dat1)
  oy_0.8 <- f.over(dat4)
  oy_0.9 <- f.over(dat7)
  
  # Optimal recruitment ----
  or_0.7 <- f.over(dat2)
  or_0.8 <- f.over(dat5)
  or_0.9 <- f.over(dat8)
  
  # Bind dataframes together
  Y <- cbind(of_0.7,oy_0.7,or_0.7,of_0.8,oy_0.8,or_0.8,of_0.9,oy_0.9,or_0.9, c(0, x))
  names(Y) <- c('of_0.7','oy_0.7','or_0.7','of_0.8','oy_0.8','or_0.8','of_0.9','oy_0.9',
                'or_0.9','Escapement')
  
  # Quantiles and Medians ----
  dat9 %>%
    summarise_all(funs(median = median, 
                       q95=quantile(., 0.95, na.rm=T), 
                       q90=quantile(., 0.90, na.rm=T),
                       q10=quantile(., 0.10, na.rm=T),
                       q5=quantile(., 0.05, na.rm=T))) -> mq
  names(mq) <- c(rep(('Median'),length(x)+1), 
                 rep(('q95'),length(x)+1), 
                 rep(('q90'),length(x)+1), 
                 rep(('q10'),length(x)+1), 
                 rep(('q5'),length(x)+1))
  
  qm <- data.frame(measure = names(mq), value = as.numeric(mq[1,]), Escapement=rep(c(0,x), length(unique(names(mq)))))
  qm <- spread(qm, measure, value)
  qm <- qm[c("q95", "q90", "Median","q10", "q5", "Escapement")]
  Y <- Y[c("oy_0.9", "oy_0.8", "oy_0.7", "of_0.9", "of_0.8", "of_0.7","or_0.9","or_0.8","or_0.7","Escapement")]
  write.csv(qm,("state_space_model/2023_analysis/update/output/base_case/processed/QM.csv"), row.names=FALSE)
  write.csv(Y,("state_space_model/2023_analysis/update/output/base_case/processed/Y.csv"), row.names=FALSE)
  
  # confidence intervals ----
  dat10 %>%
    summarise_all(funs(median = median, 
                       q95=quantile(., 0.95, na.rm=T), 
                       q90=quantile(., 0.90, na.rm=T),
                       q10=quantile(., 0.10, na.rm=T),
                       q5=quantile(., 0.05, na.rm=T))) -> mq
  names(mq) <- c(rep(('Median'),length(x)+1), 
                 rep(('q95'),length(x)+1), 
                 rep(('q90'),length(x)+1), 
                 rep(('q10'),length(x)+1), 
                 rep(('q5'),length(x)+1))
  
  CI <- data.frame(measure = names(mq), value = as.numeric(mq[1,]), Escapement=rep(c(0,x), length(unique(names(mq)))))
  CI <- spread(CI, measure, value)
  CI <- CI[c("q95", "q90", "Median","q10", "q5", "Escapement")]
  write.csv(CI,("state_space_model/2023_analysis/update/output/base_case/processed/CI.csv"), row.names=FALSE)
  read.csv("state_space_model/2023_analysis/update/output/base_case/processed/Y.csv") -> Y
  
  Y %>% 
    dplyr::select(Escapement, oy_0.9, oy_0.8, oy_0.7) %>% 
    gather(key="variable", value="value", -Escapement) %>% 
    mutate(sra = "Yield Profile",
           max_pct =ifelse(grepl("oy_0.7",variable), 
                           0.7,
                           ifelse(grepl("oy_0.8",variable),0.8, 0.9 )))-> my1
  
  Y %>% 
    dplyr::select(Escapement, of_0.9, of_0.8, of_0.7) %>% 
    gather(key="variable", value="value", -Escapement) %>% 
    mutate(sra = "Overfishing Profile",
           max_pct =ifelse(grepl("of_0.7",variable), 
                           0.7,
                           ifelse(grepl("of_0.8",variable),0.8, 0.9)))-> my2
  
  
  Y %>% 
    dplyr::select(Escapement, or_0.9, or_0.8, or_0.7) %>% 
    gather(key="variable", value="value", -Escapement) %>% 
    mutate(sra = "Recruitment Profile",
           max_pct =ifelse(grepl("or_0.7",variable), 
                           0.7,
                           ifelse(grepl("or_0.8",variable),0.8, 0.9 )))-> my3
  
  
  my4<-rbind(my1, my2, my3)
  my4 %>%
    dplyr::select(Escapement, variable, value, sra, max_pct) %>%
    mutate(Escapement = as.numeric(Escapement),
           Probability = as.numeric(value),
           max_pct = as.factor(max_pct)) -> my4
  my4 %>%
    filter(sra == "Yield Profile") -> fig_data1
  
  my4 %>%
    filter(sra == "Overfishing Profile") -> fig_data2
  
  my4 %>%
    filter(sra == "Recruitment Profile") -> fig_data3
  
  # EXPECTED SUSTAINED YIELD 
  out.file <- paste0("state_space_model/2023_analysis/update/output/base_case/processed/expect_yield.png")
  theme_set(theme_report(base_family = "Times", base_size = 12))
  ggplot(qm, aes(Escapement, Median))+geom_line(size=1)+
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha=.15)+
    geom_ribbon(aes(ymin = q10, ymax = q90), alpha=.15)+ xlab("Escapement")+
    ylab("Expected Yield")+scale_y_continuous(labels = comma)+
    scale_x_continuous(labels = comma,breaks = seq(0, 200000, 25000), limits = c(0,200000))+
    scale_y_continuous(labels = comma,breaks = seq(-200000, 200000, 25000), limits = c(-200000,200000))+
    geom_vline(xintercept = LowerB,linetype = "longdash" )+geom_vline(xintercept = UpperB ,linetype = "longdash")+
    geom_vline(xintercept = SMSY,linetype = 1 )
  ggsave(out.file, dpi = 1000, height = 4, width = 6, units = "in")
  
  # SEAK plot 
  theme_set(theme_report(base_family = "Times", base_size = 12))
  ggplot(fig_data1, aes(x = Escapement, y = Probability, linetype = max_pct)) + ggtitle("(c) Yield Profile") + 
    annotate("rect", xmin = 40000, xmax = 75000, ymin = 0, ymax = 1,
             inherit.aes = FALSE, fill = "grey80", alpha = 0.9) +
    geom_line() +    theme(plot.title = element_text(size = 12, face = "bold"),
                           strip.text.y = element_text(size=0),legend.position= "none") +
    scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
    scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1))+
    scale_linetype_discrete(name = "Percent of Max.") + xlab('Escapement (S)')+
    facet_grid(sra ~ .)  -> plot1
  
  ggplot(fig_data2, aes(x = Escapement, y = Probability, linetype = max_pct)) + 
    annotate("rect", xmin = 40000, xmax = 75000, ymin = 0, ymax = 1,
             inherit.aes = FALSE, fill = "grey80", alpha = 0.9) + ggtitle("(a) Overfishing Profile") + 
    theme(plot.title = element_text(size = 12, face = "bold"),
          strip.text.y = element_text(size=0),legend.position=c(0.92,0.86), legend.title = element_blank()) +
    geom_line() + xlab('Escapement (S)') +
    scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
    scale_linetype_discrete(name = "Percent of Max.") + 
    facet_grid(sra ~ .) -> plot2
  
  ggplot(fig_data3, aes(x = Escapement, y = Probability, linetype = max_pct)) + 
    annotate("rect", xmin = 40000, xmax = 75000, ymin = 0, ymax = 1,
             inherit.aes = FALSE, fill = "grey80", alpha = 0.9) + ggtitle("(b) Recruitment Profile") + 
    geom_line() + xlab('Escapement (S)') +   theme(plot.title = element_text(size = 12, face = "bold"),
                                                   strip.text.y = element_text(size=0),legend.position= "none") +
    scale_x_continuous(labels = comma, breaks = seq(0, 100000, 20000), limits = c(0, 100000))+
    scale_linetype_discrete(name = "Percent of Max.") +
    facet_grid(sra ~ .)  -> plot3
  cowplot::plot_grid(plot2,plot3,plot1, align = "v", nrow = 3, ncol=1) 
  ggsave("state_space_model/2023_analysis/update/output/base_case/processed/SEAK_report.png", dpi = 500, height = 7, width = 6, units = "in")}

# function to create a text file---
create_txt <- function(mu = 4, sigma = 0.075, dist = dnorm, path = "state_space_model/code/") {
  fn <- paste0("Taku_sockeye_beta_mu", mu, "_sigma", sigma, "_dist_", dist,".txt")
  txt <- c("
model
{
  for (y in (A + a.min):(Y + A - 1)) {
  log.R[y] ~ dnorm(log.R.mean2[y], tau.R)
  R[y] <- exp(log.R[y])
  log.R.mean1[y] <- log(S[y - a.max]) + lnalpha - beta * 
  S[y - a.max]
  log.resid[y] <- log(R[y]) - log.R.mean1[y]
  }
  log.R.mean2[A + a.min] <- log.R.mean1[A + a.min] + phi * 
  log.resid.0
  for (y in (A + a.min + 1):(Y + A - 1)) {
  log.R.mean2[y] <- log.R.mean1[y] + phi * log.resid[y - 
  1]
  }
  lnalpha ~ dnorm(0.00000E+00, 1.00000E-06)  T(1.00000E-06, )
  beta ~",dist,"(",mu, ",",sigma, ")
  phi ~ dnorm(0.00000E+00, 1.00000E-06)  T(-1, 1)
    mean.log.RO ~ dnorm(0.00000E+00, 1.00000E-06)
  tau.RO ~ dgamma(0.001, 0.001)
  log.resid.0 ~ dnorm(0.00000E+00, tau.red)
  tau.R ~ dgamma(0.001, 0.001)
  sigma.R <- 1/sqrt(tau.R)
  alpha <- exp(lnalpha)
  sigma.RO <- 1/sqrt(tau.RO)
  tau.red <- tau.R * (1 - phi * phi)
  lnalpha.c <- lnalpha + (sigma.R * sigma.R/2/(1 - phi * phi))
  R.O <- exp(mean.log.RO)
  for (y in 1:a.max) {
  log.R[y] ~ dnorm(mean.log.RO, tau.RO)
  R[y] <- exp(log.R[y])
  }
  S.max <- 1/beta
  alpha.c <- min(exp(lnalpha.c), 10000)
  S.eq.c <- lnalpha.c * S.max
  U.msy.c <- lnalpha.c * (0.5 - 0.07 * lnalpha.c)
  S.msy.c <- S.eq.c * (0.5 - 0.07 * lnalpha.c)
  positive.lna.c <- step(lnalpha.c)
  lnalpha.c.nonneg <- lnalpha.c * positive.lna.c
  S.eq.c2 <- lnalpha.c.nonneg * S.max
  peterman.approx.c <- (0.5 - 0.65 * pow(lnalpha.c.nonneg, 
  1.27)/(8.7 + pow(lnalpha.c.nonneg, 1.27)))
  U.msy.c2 <- lnalpha.c.nonneg * peterman.approx.c
  S.msy.c2 <- U.msy.c2/beta
  U.max.c2 <- 1 - 1/exp(lnalpha.c.nonneg)
  S.msy.c.80 <- S.msy.c * 0.8
  S.msy.c2.80 <- S.msy.c2 * 0.8
  D.scale ~ dunif(0.00000E+00, 1)
  D.sum <- 1/(D.scale * D.scale)
  pi.2p ~ dbeta(1, 1)
  pi.1 ~ dbeta(1, 1)
  pi[1] <- pi.1
  pi[2] <- pi.2p * (1 - pi[1])
  pi[3] <- 1 - pi[1] - pi[2]
  for (a in 1:A) {
  gamma[a] <- D.sum * pi[a]
  for (y in 1:(Y + A - 1)) {
  g[y, a] ~ dgamma(gamma[a], 0.01)
  p[y, a] <- g[y, a]/sum(g[y, ])
  }
  }
  for (a in 1:A) {
  for (y in a:(Y + (a - 1))) {
  N.ya[y - (a - 1), (A + 1 - a)] <- p[y, (A + 1 - a)] * 
  R[y]
  }
  }
  for (y in 1:Y) {
  N[y] <- sum(N.ya[y, 1:A])
  for (a in 1:A) {
  q[y, a] <- N.ya[y, a]/N[y]
  }
  }
  for (t in 1:Y) {
  x[t, 1:A] ~ dmulti(q[t, ], n.a[t])
  }
  for (y in 1:Y) {
  mu.hbelow_wild[y] ~ dbeta(1, 1)
  h.below_wild[y] <- mu.hbelow_wild[y] * N[y]
  log.hb_wild[y] <- log(h.below_wild[y])
  tau.log.hb.wild[y] <- 1/log(cv.hb[y] * cv.hb[y] + 1)
  hbelow_wild[y] ~ dlnorm(log.hb_wild[y], tau.log.hb.wild[y])
  mu.habove_wild[y] ~ dbeta(1, 1)
  h.above_wild[y] <- mu.habove_wild[y] * N[y]
  log.ha_wild[y] <- log(h.above_wild[y])
  tau.log.ha.wild[y] <- 1/log(cv.ha[y] * cv.ha[y] + 1)
  habove_wild[y] ~ dlnorm(log.ha_wild[y], tau.log.ha.wild[y])
  inriver.run[y] <- max(N[y] - h.below_wild[y], 1)
  log.ir[y] <- log(inriver.run[y])
  tau.log.ir[y] <- 1/log(cv.ir[y] * cv.ir[y] + 1)
  ir[y] ~ dlnorm(log.ir[y], tau.log.ir[y])
  mu.habove[y] ~ dbeta(1, 1)
  h.above[y] <- mu.habove[y] * inriver.run[y]
  log.ha[y] <- log(h.above[y])
  tau.log.ha[y] <- 1/log(cv.ha[y] * cv.ha[y] + 1)
  habove[y] ~ dlnorm(log.ha[y], tau.log.ha[y])
  mu.hbelow[y] ~ dbeta(1, 1)
  h.below[y] <- mu.hbelow[y] * N[y]
  log.hb[y] <- log(h.below[y])
  tau.log.hb[y] <- 1/log(cv.hb[y] * cv.hb[y] + 1)
  hbelow[y] ~ dlnorm(log.hb[y], tau.log.hb[y])
  mu[y] <- (h.below[y] + h.above[y])/N[y]
  S[y] <- max(inriver.run[y] - h.above[y], 1)
  log.S[y] <- log(S[y])
  }
}")
  
  write.table(txt, file = paste0(path, fn),  quote = FALSE, sep=" ",
              row.names = FALSE, col.names = FALSE)
}

