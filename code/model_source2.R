#State Space Model Function with autocorrelation
jag_model_SR <- function(){
  
  for (y in 1:nyrs){
    # log normal Likelihood 
    R[y] ~ dlnorm(mu[y],Tau)
    # rk = 1 and bh = 0:  log-linearlized Rikcer SR Model 
    #      fit[y] <- log(S[y]) + lnalpha - rk*(beta * S[y]/(10^d))
    # rk = 0 and bh = 1:  log-linearlized Beverton-Holt Model 
    #      fit[y] <- log(S[y]) + lnalpha - (log(1+beta*S[y]/(10^d)))   
    fit[y] <- log(S[y]) + lnalpha - rk*(beta * S[y]/(10^d))-bh*(log(1+beta*S[y]/(10^d)))
    e[y] <- log(R[y]) - fit[y]
  }
  # ar1 = 0 and kf =0 in standard analysis   
  # ar1 = 1 and kf =0 AR1 error moddel 
  # ar1 = 0 and kf =1 time-varying alpha
  # e0 original residual     
  mu[1] <-  fit[1] + ar1*phi * e0;
  cw[1]  ~ dnorm(0,tauw)
  for (y in 2:nyrs){	   
    cw[y] ~  dnorm(cw[y-1],tauw)
    mu[y] <- fit[y] + kf*cw[y] + ar1*phi*e[y-1]
  }
  #   Define Priors
  lnalpha ~ dunif(min.a,max.a)
  beta ~ dunif(min.b,max.b)
  sigma ~ dunif(0,2)  
  sigmaw ~ dunif(0,2)
  #    sigmav ~ dunif(0,10)
  phi ~ dnorm(0,5)%_%T(-1,1) 
  e0 ~ dnorm(0,25) 
  Tau <- 1/(sigma*sigma)
  tauw <- 1/(sigmaw*sigmaw)	
  #   tauv <- 1/(sigmav*sigmav)
  # Extract time-varyihg alapha
  for(y in 1:nyrs){
    lnalphai[y] <- lnalpha+cw[y] 
  }
}