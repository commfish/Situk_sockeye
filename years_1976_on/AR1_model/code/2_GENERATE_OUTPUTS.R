# THIS SCRIPT IS SOURCED FROM INSIDE 1_RUN_MODEL.R
#(if package.use = "rjags") 

#Gelman statistic
#Brooks and Gelman (1997) have suggested, if Rc<1.2 for all model parameters, 
#one can be fairly confident that convergence has been reached. Otherwise, longer chains or other means for improving the convergence may be needed
#Brooks, S. P., and A. Gelman. 1997. General Methods for Monitoring Convergence of Iterative Simulations. 
#Journal of Computational and Graphical Statistics 7: 434–455.
d <- 4
windowsFonts(Times=windowsFont("Arial"))
theme_set(theme_sleek())
out.path <- paste0("years_1976_on/AR1_model/output")
parameters <- c("lnalpha", "phi", "beta", "sigma", "sigmaw", "Tau", "tauw", "alpha", "lnalpha.c", "e0")

# Numerical summary of each parameter (mean, median, quantiles of posteriors)----
# https://stats.stackexchange.com/questions/429470/what-is-the-correct-effective-sample-size-ess-calculation
coda::effectiveSize(post) %>%
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
g <- geweke.diag(post)

zmat <- do.call(cbind, lapply(g, function(x) x$z))# z=0 is good convergence; >2 is potential issue

geweke_df <- as.data.frame(zmat)
geweke_df$variable <- rownames(geweke_df)
rownames(geweke_df) <- NULL

geweke_df <- dplyr::relocate(geweke_df, variable)
colnames(geweke_df)[2:ncol(geweke_df)] <- paste0("chain", seq_len(ncol(geweke_df) - 1))
write.csv(geweke_df, paste0(out.path,"/geweke.csv")) 

pdf("years_1976_on/AR1_model/output/geweke.pdf",height=10, width=8,onefile=T)
geweke.plot(post)
dev.off()

# 90th and 95th percentiles of output quants----
get_quant <- function(parameters) {
  x <- as.vector(post.arr[, parameters, ])   # extract all chains
  q <- quantile(x, probs = c(0, .025, .05, .5, .95, .975, 1))
  data.frame(parameter = parameters, t(q))
}

quant_table <- do.call(rbind, lapply(parameters, get_quant))

write.csv(quant_table, paste0(out.path, "/percentiles.csv"), row.names = FALSE)

# lambert calc----
# Extract parameters from posterior array
parameters <- c("lnalpha", "lnalpha.c","beta")

# Combine chains into single data frame
coda <- as.data.frame(post.arr[, parameters, ]) %>%
  {
    coda1 <- .[, 1:3]; coda2 <- .[, 4:6]; coda3 <- .[, 7:9]
    names(coda1) <- parameters
    names(coda2) <- parameters
    names(coda3) <- parameters
    rbind(coda1, coda2, coda3)
  } %>%
  mutate(
    beta1 = beta * 10^-4,
    
# alpha
    alpha = exp(lnalpha),
    alpha.c   = exp(lnalpha.c),

# Smsy - spawners at maximum sustainable yield
Smsy_lambert   = (10^d) * (1 - lambert_W0(exp(1 - lnalpha)))   / beta,
Smsy_lambert.c = (10^d) * (1 - lambert_W0(exp(1 - lnalpha.c))) / beta,

# Umsy - exploitation rate at MSY
Umsy_lambert   = Smsy_lambert   * beta / (10^d),
Umsy_lambert.c = Smsy_lambert.c * beta / (10^d),  

# Smax - spawners at maximum recruitment
Smax = (10^d) / beta,
    
# Rmsy - recruits at MSY
Rmsy.c = Smsy_lambert.c * exp(lnalpha.c - beta1 * Smsy_lambert.c),
Rmsy   = Smsy_lambert   * exp(lnalpha   - beta1 * Smsy_lambert),
 
# Seq - equilibrium spawner abundance
Seq   = lnalpha / beta1,

# MSY - maximum sustainable yield
MSY   = Rmsy   - Smsy_lambert,
MSY.c   = Rmsy.c   - Smsy_lambert.c,    

# Rmax - maximum recruitment
Rmax.c = exp(lnalpha.c) * (1 / beta1) * exp(-1),
Rmax   = exp(lnalpha)   * (1 / beta1) * exp(-1)) %>%
  as.data.frame()

write.csv(coda, file= paste0(out.path,"/coda.csv") ,row.names=FALSE)  

qtab <- apply(coda, 2, quantile, #quantiles
              probs = c(0, 0.025, 0.05, 0.5, 0.95, 0.975, 1)) %>%
  as.data.frame() %>%
  rownames_to_column("stat")

mtab <- apply(coda, 2, mean) %>% #means
  t() %>%
  as.data.frame() %>%
  mutate(stat = "mean")

stab <- apply(coda, 2, sd) %>% # stdev
  t() %>%
  as.data.frame() %>%
  mutate(stat = "sd")

x <- bind_rows(qtab, mtab, stab) %>% # combine
  column_to_rownames("stat") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("variable")

names(x) <- c("variable","0","2.5","5","50","95","97.5","100","mean","sd")

x %>%
  mutate(across(-variable, ~as.numeric(.))) %>%
  mutate(post_cv = sd/mean) %>%
  write.csv(., file= paste0(out.path,"/quantiles_lambert.csv"))    

# lambert density plot setup----
post_mat <- as.matrix(post)
lnalpha_draws  <- post_mat[, "lnalpha"]
lnalphac_draws <- post_mat[, "lnalpha.c"]
beta_draws     <- post_mat[, "beta"]

Smsy   <- (1 - lambert_W0(exp(1 - lnalpha_draws))) / (beta_draws*10^-4)
Smsy.c <- (1 - lambert_W0(exp(1 - lnalphac_draws))) / (beta_draws*10^-4)
Umsy   <-  1 - lambert_W0(exp(1 - lnalpha_draws))


df1 <- data.frame(Smsy = Smsy)
df2 <- data.frame(Umsy = Umsy)

# lambert density plots----
out.file <- paste0("years_1976_on/AR1_model/output/density.png")
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
  mutate(iter = 1) %>%
  pivot_longer(
    cols = -iter,
    names_to = "parameter",
    values_to = "value")
out.file <- paste0("years_1976_on/AR1_model/output/posterior_densities.png")
ggplot(post_df, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.4, color = "black") +
  facet_wrap(~ parameter, scales = "free", ncol = 3) +
  theme_bw() +
  labs(x = "Value", y = "Density", title = "Posterior densities")
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

# trace plots
df_trace <- as.data.frame.table(post.arr, responseName = "value") %>%
  dplyr::rename(
    iter  = iter,
    var = var,
    chain = chain,
    value  = value) %>%
  mutate(
    iter  = as.numeric(iter),
    chain = as.factor(chain))

out.file <- paste0("years_1976_on/AR1_model/output/trace_plots.png")
ggplot(df_trace, aes(x = iter, y = value, color = chain)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  facet_wrap(~ var, scales = "free_y") +
  labs(x = "Iteration", y = "Value", color = "Chain") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom")
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

out.file <- paste0("years_1976_on/AR1_model/output/dens_plots.png")
dens_plot <- ggplot(df_trace, aes(value, fill = chain)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ var, scales = "free") +
  theme_bw()
ggsave(out.file, dpi = 500, height = 8, width = 9, units = "in")

# recruit residuals
parameters <- c("lnalpha", "beta", "lnalpha.c", "e0", "phi")
x <- as.data.frame(post.arr[, parameters, ])

coda1 <- x[, 1:5]
coda2 <- x[, 6:10]
coda3 <- x[, 11:15]

names(coda1) <- parameters
names(coda2) <- parameters
names(coda3) <- parameters

rbind(coda1, coda2, coda3) %>%
  as.data.frame() -> SR.post

# recruit plot data
read.csv("data/Situk_sockeye.csv") %>%
  mutate(yield = (recruit50 - spawn),
         logR = log(recruit50)) -> dat
R <- dat$recruit50
S <- dat$spawn
d <- floor(log10(mean(S)))

# MCMC parameters
lnalpha <-SR.post$lnalpha
beta <- SR.post$beta
phi <- SR.post$phi
e0 <- SR.post$e0

# dimensions
n_years <- length(S)
n_mcmc  <- length(lnalpha)

# empty matrices
RD  <- matrix(NA, nrow = n_mcmc, ncol = n_years)  # basic Ricker residuals
RD2 <- matrix(NA, nrow = n_mcmc, ncol = n_years)  # AR1 corrected residuals

# basic Ricker residuals
for(i in 1:n_mcmc){
  Ey      <- srmodel(lnalpha[i], beta[i], S, d)
  RD[i,]  <- log(R) - log(Ey)
}

# AR1 corrected residuals
for(i in 1:n_mcmc){
  RD2[i, 2:n_years] <- RD[i, 2:n_years] - phi[i] * RD[i, 1:(n_years - 1)]
  RD2[i, 1] <- RD[i, 1] - phi[i] * e0[i]
}

# convert to data frames with named columns
RD_df  <- as.data.frame(RD)
RD2_df <- as.data.frame(RD2)

names(RD_df)  <- paste0("RD.", 1:n_years)
names(RD2_df) <- paste0("RD2.", 1:n_years)

# Combine and write to CSV
out <- cbind(RD_df, RD2_df)
write.csv(out, file = "years_1976_on/AR1_model/output/processed/recruit_data.csv", row.names = FALSE)

# create param table
read.csv("years_1976_on/AR1_model/output/quantiles_lambert.csv") %>%
  dplyr::select(variable, X2.5, X50, X97.5, post_cv) %>%
  write.csv(., file = "years_1976_on/AR1_model/output/processed/params.csv", row.names = FALSE)

# ACF plots
pdf("years_1976_on/AR1_model/output/ACF_parameters.pdf", width = 7, height = 5)
for(p in parameters){
  acf_param(post, p, combine = TRUE)
}

dev.off()

# ACF values
post_mat <- as.matrix(post)

e_cols <- grep("e[", colnames(post_mat), fixed = TRUE, value = TRUE)
e_mean <- colMeans(post_mat[, e_cols])

acf_e <- acf(e_mean, plot = FALSE)

data.frame(lag = as.numeric(acf_e$lag),acf = as.numeric(acf_e$acf))
