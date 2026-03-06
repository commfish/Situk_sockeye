# THIS SCRIPT IS SOURCED FROM INSIDE 1_RUN_MODEL.R
#(if package.use = "rjags") 
# THIS SCRIPT IS SOURCED FROM INSIDE 1_RUN_MODEL.R
#(if package.use = "rjags" and location of file is rjags_Explore_BaseCase) 

# - post.samp and post 
# these are the mcmc.list objects created by coda.samples
# you can access individual variables like this: post[,"var.name"]
# - post.arr 
# - coda
# create coda file
#Gelman statistic
#Brooks and Gelman (1997) have suggested, if Rc<1.2 for all model parameters, 
#one can be fairly confident that convergence has been reached. Otherwise, longer chains or other means for improving the convergence may be needed
#Brooks, S. P., and A. Gelman. 1997. General Methods for Monitoring Convergence of Iterative Simulations. 
#Journal of Computational and Graphical Statistics 7: 434–455.
#loadfonts(device="win")
windowsFonts(Times=windowsFont("Arial"))
theme_set(theme_sleek())
out.path <- paste0("output/", out.label)
parameters <- c("lnalpha", "phi", "beta", "sigma", "sigmaw", "Tau", "tauw", "alpha", "lnalpha.c")

# Numerical summary of each parameter (mean, median, quantiles of posteriors)----
effectiveN<-coda::effectiveSize(post) #https://stats.stackexchange.com/questions/429470/what-is-the-correct-effective-sample-size-ess-calculation
as.data.frame(effectiveN)%>%
  write.csv(., file= paste0(out.path,"/effectiveN.csv")) 
summary<-summary(post)  
stats<-summary$statistics;  colnames(stats)
quants<-summary$quantiles;  colnames(quants)
statsquants <- cbind(stats,quants) 
statsquants %>% 
  as.data.frame() %>% 
  dplyr::select(Mean, SD, 'Time-series SE', '2.5%', '50%', '97.5%') %>%
  dplyr::rename(time_series_se = 'Time-series SE') %>%
  rownames_to_column('variable') %>%
  mutate (mc_error = time_series_se/SD,
          converge = ifelse(mc_error < 0.05, "true", "false")) %>%
  write.csv(., file= paste0(out.path,"/statsquants.csv"))    

# Gelman Diagnostics----
gel <- as.data.frame(gelman.diag(post, multivariate=F)[[1]])
poor.threshold = 1.01 #values less than 1.01 are generally considered converged
gel %>%
  rownames_to_column('variable') %>%
  rename(point_estimate = "Point est.") %>%
  mutate (converge = ifelse(point_estimate < poor.threshold, "true", "false")) %>%
  write.csv(., file= paste0(out.path,"/gelman.csv"))   

# Geweke Diagnostics----
#Examine convergence of the Markov chains using the Geweke's convergence diagnostic
#a convergence diagnostic for Markov chains based on a test for equality 
#of the means of the ﬁrst and last part of a Markov chain (by default the ﬁrst 10% and the last 50%).
#As for a cut-off you can compare to the standard normal critical values z α/2 where α=0.05
#null hypothesis is rejected if Z large (used to determine the burn-in period)
# Geweke diagnostic
l <- geweke.diag(post)
df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T),stringsAsFactors=FALSE)
df <- t(df) %>%
  as.data.frame() %>%
  rownames_to_column('variable1') -> df
colnames(df) <- c("variable1", "chain1", "chain2", "chain3")
names <-read.csv(file = paste0(out.path,"/statsquants.csv"))
names  %>%
  dplyr::select(variable) %>% 
  add_row(variable = "frac1") %>% #0.1 and 0.5 default
  add_row(variable = "frac2")  -> names
x <- cbind(names, df)
x  %>% 
  dplyr::select(-variable1) %>%
  write.csv(., paste0(out.path,"/geweke.csv")) 

pdf("output/base_case/geweke.pdf",height=10, width=8,onefile=T)
geweke.plot(post)
dev.off()

# 90th and 95th percentiles of output quants----
parameters=c('lnalpha')
x <- post.arr[,parameters,]
lnalpha<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
lnalpha <- data.frame(lnalpha)

parameters=c('lnalpha.c')
x <- post.arr[,parameters,]
lnalpha.c<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
lnalpha.c <- data.frame(lnalpha.c)

parameters=c('alpha')
x <- post.arr[,parameters,]
alpha<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
alpha <- data.frame(alpha)

parameters=c('beta')
x <- post.arr[,parameters,]
beta<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
beta <- data.frame(beta)

parameters=c('phi')
x <- post.arr[,parameters,]
phi<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
phi <- data.frame(phi)

parameters=c('sigma')
x <- post.arr[,parameters,]
sigma<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
sigma <- data.frame(sigma)

parameters=c('sigmaw')
x <- post.arr[,parameters,]
sigmaw<-quantile(x, probs=c(0,0.025, 0.05, 0.5, 0.95, 0.975, 1))
sigmaw <- data.frame(sigmaw)

parameters=c('Tau')
x <- post.arr[,parameters,]
Tau<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
Tau <- data.frame(Tau)

parameters=c('tauw')
x <- post.arr[,parameters,]
tauw<-quantile(x, probs=c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1))
tauw<- data.frame(tauw)

step1 <- cbind(alpha, lnalpha)
step2 <- cbind(step1, lnalpha.c)
step3 <- cbind(step2, alpha)
step4 <- cbind(step3, beta)
step5 <- cbind(step4, phi)
step6 <- cbind(step5, sigma)
step7 <- cbind(step6, sigmaw)
step8 <- cbind(step7, Tau)
step9 <- cbind(step8, tauw)
step9 %>% 
  t()%>%
  write.csv(., file= paste0(out.path,"/percentiles.csv")) 


# lambert calc----
parameters=c("lnalpha", "beta", "lnalpha.c")
x <- as.data.frame(post.arr[,parameters,])
coda1 <- x[,1:3]
coda2 <- x[,4:6]
coda3 <- x[,7:9]
coda1 %>% 
  dplyr::rename(beta = beta.1,
                lnalpha = lnalpha.1,
                lnalpha.c = lnalpha.c.1)-> coda1
coda2 %>% 
  dplyr::rename(beta = beta.2,
                lnalpha = lnalpha.2,
                lnalpha.c = lnalpha.c.2) -> coda2
coda3 %>% 
  dplyr::rename(beta = beta.3,
                lnalpha = lnalpha.3,
                lnalpha.c = lnalpha.c.3) -> coda3
coda<-rbind(coda1,coda2,coda3)
write.csv(coda, file= paste0(out.path,"/coda.csv") ,row.names=FALSE)  
coda %>% 
  mutate(Smsy_lambert.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, # mean recruitment
         Smsy_lambert = (1-lambert_W0(exp(1-lnalpha)))/beta, # median recruitment
         Umsy_lambert = (1-lambert_W0(exp(1-lnalpha.c))),
         Rmsy = Smsy_lambert*exp(lnalpha-beta*Smsy_lambert),
         MSY = Rmsy-Smsy_lambert)  %>%
  as.data.frame() %>%
  dplyr::select(Smsy_lambert.c,Smsy_lambert, Umsy_lambert, MSY) -> coda
coda %>% 
  apply(., 2, sd) %>%
  as.data.frame()%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('variable') %>%
  mutate(variable = ifelse(variable == '.', "sd", "sd"))-> sd
coda %>% 
  apply(., 2, mean) %>%
  as.data.frame() %>%
  t()%>%
  as.data.frame() %>%
  rownames_to_column('variable') %>%
  mutate(variable = ifelse(variable == '.', "mean", "mean"))-> mean
rbind(sd,mean) -> x1
q1<-apply(coda,2,quantile,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
q1 %>%
  as.data.frame() %>%
  rownames_to_column('variable') -> x2
rbind(x1,x2) %>%
  t() %>%
  write.csv(., file= paste0(out.path,"/quantiles_lambert.csv"))    

# lambert density plot setup----
parameters=c("lnalpha", "beta", "lnalpha.c")
x <- as.data.frame(post.arr[,parameters,])
coda1 <- x[,1:3]
coda2 <- x[,4:6]
coda3 <- x[,7:9]
coda1 %>% 
  dplyr::rename(beta = beta.1,
                lnalpha = lnalpha.1,
                lnalpha.c = lnalpha.c.1)-> coda1
coda2 %>% 
  dplyr::rename(beta = beta.2,
                lnalpha = lnalpha.2,
                lnalpha.c = lnalpha.c.2) -> coda2
coda3 %>% 
  dplyr::rename(beta = beta.3,
                lnalpha = lnalpha.3,
                lnalpha.c = lnalpha.c.3) -> coda3
coda<-rbind(coda1,coda2,coda3)
coda %>% 
  mutate(Smsy = (1-lambert_W0(exp(1-lnalpha.c)))/beta,
         Umsy = (1-lambert_W0(exp(1-lnalpha.c))))%>%
  as.data.frame() -> coda
coda %>% 
  dplyr::select(Smsy) -> Smsy
coda %>% 
  dplyr::select(Umsy) -> Umsy

# lambert density plot----
out.file <- paste0("output/base_case/density.png")
options(scipen=999)
ggplot(Smsy, aes(x=Smsy, fill=Smsy, color = Smsy)) +
  geom_density(fill ="#999999", alpha=0.5)+ annotate("text",x = 0, y=0.00004, label="a)", family="Arial" ,size=6) +
  scale_color_manual(values=c("#999999"))+
  scale_fill_manual(values=c("#999999"))+
  geom_vline(xintercept = 43857,linetype = "longdash" ) +
  labs(x="Smsy",y="Density") + theme_set(theme_bw(base_size=14,base_family=
                                                    'Arial')+
                                           theme(panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank())) + theme(legend.position="none")+ 
  scale_x_continuous(labels = comma,breaks = seq(0, 150000, 25000), limits = c(0, 150000)) -> plot1


ggplot(Umsy, aes(x=Umsy, fill=Umsy)) +
  geom_density(fill ="#999999", alpha=0.5)+ annotate("text",x = 0, y=5, label="c)", family="Arial" ,size=6) +
  scale_color_manual(values=c("#999999"))+
  scale_fill_manual(values=c("#999999"))+geom_vline(xintercept = 0.75,linetype = "longdash" ) +
  labs(x="Umsy",y="Density") + theme_set(theme_bw(base_size=14,base_family=
                                                    'Arial')+
                                           theme(panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank())) + theme(legend.position="none")+ 
  scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) -> plot2
cowplot::plot_grid(plot1,plot2,  align = "v", nrow = 3, ncol=1) 
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

#trace and density plots----
# parameters <- c("lnalpha", "phi", "beta", "sigma", "sigmaw", "Tau", "tauw", "alpha")
# png("output/base_case/density_other.png",res=500, width=8, height=9, units ="in")
# denplot(post, parms = c(parameters), style="plain", greek=TRUE, col="gray30", collapse =T)
# dev.off()
# dev.off()
# 
# parameters <- c('S')
# pdf("output/base_case/density2.pdf",height=10, width=8,onefile=F)
# denplot(post, parms = c(parameters))
# dev.off()
# pdf("output/base_case/trace2.pdf",height=10, width=8,onefile=F)
# traplot(post, parms = c(parameters))
# dev.off()
# parameters <- c('R')
# pdf("output/base_case/density3.pdf",height=10, width=8,onefile=F)
# denplot(post, parms = c(parameters))
# dev.off()
# pdf("output/base_case/trace3.pdf",height=10, width=8,onefile=F)
# traplot(post, parms = c(parameters))
# dev.off()

# autocorrelation plots----
windows(record=T)
pdf("output/base_case/autocorr.pdf",height=6, width=8,onefile=T,useDingbats=F)
autocorr.plot(post, lag.max=5)
dev.off()
dev.off()
autocorr.summary<-autocorr.diag(post)
autocorr.summary<-data.frame(autocorr.summary)
write.csv(autocorr.summary, file= paste0(out.path,"/autocorr.csv")) 


# #density and time series plots----
# post.samp <- post
# nvars  <- dim(post.samp[[1]])[2]
# nsamps <- dim(post.samp[[1]])[1]
# int    <- 30
# 
# pdf("output/base_case/profiles.pdf", height=6, width=8)
# 
# for (j in seq(1, nvars, int)) {
#   
#   par(mfrow=c(5,4), mai=c(0.3,0.3,0.2,0.2))
#   
#   for (i in 0:(int-1)) {
#     
#     idx <- i + j
#     if (idx > nvars) break
#     
#     # Combine all chains for consistent axis limits
#     allvals <- c(post.samp[[1]][,idx],
#                  post.samp[[2]][,idx],
#                  post.samp[[3]][,idx])
#     
#     mindat <- min(allvals)
#     maxdat <- max(allvals)
#     
#     # --- Density plot ---
#     plot(density(post.samp[[1]][,idx]), col='blue',
#          main=colnames(post.samp[[1]])[idx],
#          xlim=c(mindat, maxdat))
#     lines(density(post.samp[[2]][,idx]), col='red')
#     lines(density(post.samp[[3]][,idx]), col='green')
#     
#     # --- Trace plot ---
#     plot(post.samp[[1]][,idx], type='l', col='blue',
#          main=colnames(post.samp[[1]])[idx],
#          ylim=c(mindat, maxdat))
#     lines(post.samp[[2]][,idx], col='red')
#     lines(post.samp[[3]][,idx], col='green')
#   }
# }
# 
# dev.off()
# 
# 
