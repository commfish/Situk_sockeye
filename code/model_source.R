#State Space Model Function with autocorrelation
jag.model.SR = function(){
  for(y in 1:nyrs){
    # log normal Likelihood 
    R[y] ~ dlnorm(mu[y],Tau)
    fit[y] <- log(S[y]) + lnalpha - (beta * S[y])
    e[y] <- log(R[y]) - fit[y]
  }
    
  mu[1] <-  fit[1] + phi*e0
  for(y in 2:nyrs){	   
    mu[y] <- fit[y] + phi*e[y-1]
  }
  #   Define Priors
  lnalpha ~ dunif(0,10)
  beta ~ dunif(0,10)
  sigma ~ dunif(0,2)  
  sigmaw ~ dunif(0,2)
  phi ~ dnorm(0,5)%_%T(-1,1) 
  e0 ~ dnorm(0,25) 
  Tau <- 1/(sigma*sigma)
  tauw <- 1/(sigmaw*sigmaw)	
  alpha <- exp(lnalpha)
  lnalpha.c <-lnalpha+0.5*(sigma^2)/(1-phi^2) # Calculate Mean R adjustment    
  alpha.c <-exp(lnalpha.c)
  }

