# THIS SCRIPT IS SOURCED FROM INSIDE 1_RUN_MODEL.R
#(if package.use = "rjags") 

#Gelman statistic
#Brooks and Gelman (1997) have suggested, if Rc<1.2 for all model parameters, 
#one can be fairly confident that convergence has been reached. Otherwise, longer chains or other means for improving the convergence may be needed
#Brooks, S. P., and A. Gelman. 1997. General Methods for Monitoring Convergence of Iterative Simulations. 
#Journal of Computational and Graphical Statistics 7: 434–455.
windowsFonts(Times=windowsFont("Arial"))
theme_set(theme_sleek())
out.path <- paste0("output/", out.label)
parameters <- c("lnalpha", "phi", "beta", "sigma", "sigmaw", "Tau", "tauw", "alpha", "lnalpha.c")

# Numerical summary of each parameter (mean, median, quantiles of posteriors)----
coda::effectiveSize(post) %>%#https://stats.stackexchange.com/questions/429470/what-is-the-correct-effective-sample-size-ess-calculation
  as.data.frame %>%
  write.csv(., file= paste0(out.path,"/effectiveN.csv")) 

summary <- summary(post)  
stats<-summary$statistics;  colnames(stats)
quants<-summary$quantiles;  colnames(quants)
cbind(stats,quants) %>% 
  as.data.frame() %>% 
  dplyr::select(Mean, SD, 'Time-series SE', '2.5%', '50%', '97.5%') %>%
  dplyr::rename(time_series_se = 'Time-series SE') %>%
  rownames_to_column('variable') %>%
  mutate (mc_error = time_series_se/SD,
          converge = ifelse(mc_error < 0.05, "true", "false")) %>%
  write.csv(., file= paste0(out.path,"/statsquants.csv"))    

# Gelman Diagnostics----
poor.threshold <- 1.01  # R-hat < 1.01 = converged

gel_df <- gelman.diag(post, multivariate = FALSE)[[1]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  dplyr::rename(point_estimate = `Point est.`) %>%
  dplyr::mutate(converge = ifelse(point_estimate < poor.threshold, "true", "false"))

write.csv(gel_df, file = file.path(out.path, "gelman.csv"), row.names = FALSE)

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

# compute Geweke diagnostics
g <- geweke.diag(post)

zmat <- do.call(cbind, lapply(g, function(x) x$z))# z=0 is good convergence; >2 is potential issue

geweke_df <- as.data.frame(zmat)
geweke_df$variable <- rownames(geweke_df)
rownames(geweke_df) <- NULL

geweke_df <- dplyr::relocate(geweke_df, variable)
colnames(geweke_df)[2:ncol(geweke_df)] <- paste0("chain", seq_len(ncol(geweke_df) - 1))
geweke_df

pdf("output/base_case/geweke.pdf",height=10, width=8,onefile=T)
geweke.plot(post)
dev.off()

# 90th and 95th percentiles of output quants----
params <- c("lnalpha", "lnalpha.c", "alpha", "beta", "phi", "sigma", "sigmaw", "Tau", "tauw")

get_quant <- function(param) {
  x <- as.vector(post.arr[, param, ])   # extract all chains
  q <- quantile(x, probs = c(0, .025, .05, .5, .95, .975, 1))
  data.frame(parameter = param, t(q))
}

quant_table <- do.call(rbind, lapply(params, get_quant))

write.csv(quant_table, paste0(out.path, "/percentiles.csv"), row.names = FALSE)

# lambert calc----
parameters <- c("lnalpha", "beta", "lnalpha.c")
x <- as.data.frame(post.arr[, parameters, ])

coda1 <- x[, 1:3]
coda2 <- x[, 4:6]
coda3 <- x[, 7:9]

names(coda1) <- parameters
names(coda2) <- parameters
names(coda3) <- parameters

coda <- rbind(coda1, coda2, coda3)
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
  mutate(variable = ifelse(variable == '.', "mean", "mean")) %>%
rbind(., sd) -> x1

q1 <- apply(coda,2,quantile,probs=c(0,0.025,0.05,0.5,0.95,0.975,1))
q1 %>%
  as.data.frame() %>%
  rownames_to_column('variable') -> x2
rbind(x1,x2) %>%
  t() %>%
  write.csv(., file= paste0(out.path,"/quantiles_lambert.csv"))    

# lambert density plot setup----
post_mat <- as.matrix(post)
lnalpha_draws  <- post_mat[, "lnalpha"]
lnalphac_draws <- post_mat[, "lnalpha.c"]
beta_draws     <- post_mat[, "beta"]

Smsy   <- (1 - lambert_W0(exp(1 - lnalpha_draws))) / beta_draws
Smsy.c <- (1 - lambert_W0(exp(1 - lnalphac_draws))) / beta_draws
Umsy   <-  1 - lambert_W0(exp(1 - lnalpha_draws))


df1 <- data.frame(Smsy = Smsy)
df2 <- data.frame(Umsy = Umsy)

# lambert density plots----
out.file <- paste0("output/base_case/density.png")
theme_set(theme_bw(base_size = 14, base_family = "Arial"))

plot1 <- ggplot(df1, aes(x = Smsy)) +
  geom_density(fill = "#999999", color = "#999999", alpha = 0.5) +
  annotate("text", x = 0, y = 0.00004, label = "a)", family = "Arial", size = 6) +
  geom_vline(xintercept = 43857, linetype = "longdash") +
  labs(x = "Smsy", y = "") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none") +
  scale_x_continuous(
    labels = scales::comma,
    breaks = seq(0, 150000, 25000),
    limits = c(0, 150000))


plot2 <- ggplot(df2, aes(x = Umsy)) +
  geom_density(fill = "#999999", color = "#999999", alpha = 0.5) +
  annotate("text", x = 0, y = 5, label = "b)", family = "Arial", size = 6) +
  geom_vline(xintercept = 0.75, linetype = "longdash") +
  labs(x = "Umsy", y = "") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none") +
  scale_x_continuous(
    breaks = seq(0, 1, 0.25),
    limits = c(0, 1))

plot1 <- plot1 +
  theme(
    axis.title.y = element_text(margin = margin(r = 12)),
    plot.margin = margin(l = 20, r = 10, t = 10, b = 10))

plot2 <- plot2 +
  theme(
    axis.title.y = element_text(margin = margin(r = 12)),
    plot.margin = margin(l = 20, r = 10, t = 10, b = 10))

cowplot::plot_grid(plot1,plot2,  align = "v", nrow = 2, ncol=1) 
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

# post is your mcmc.list from coda.samples()
post_df <- as.data.frame(as.mcmc(do.call(rbind, post))) %>%
  mutate(iter = row_number()) %>%
  pivot_longer(
    cols = -iter,
    names_to = "parameter",
    values_to = "value")
out.file <- paste0("output/base_case/posterior_densities.png")
ggplot(post_df, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.4, color = "black") +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  theme_bw() +
  labs(x = "Value", y = "Density", title = "Posterior densities")
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

# trace plots
df_trace <- as.data.frame.table(post.arr, responseName = "value") %>%
  rename(
    iter  = iter,
    param = var,
    chain = chain) %>%
  mutate(
    iter  = as.numeric(iter),
    chain = as.factor(chain))

out.file <- paste0("output/base_case/trace_plots.png")
ggplot(df_trace, aes(x = iter, y = value, color = chain)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  facet_wrap(~ param, scales = "free_y") +
  labs(x = "Iteration", y = "Value", color = "Chain") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom")
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

out.file <- paste0("output/base_case/dens_plots.png")
dens_plot <- ggplot(df_trace, aes(value, fill = chain)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ param, scales = "free") +
  theme_bw()
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

