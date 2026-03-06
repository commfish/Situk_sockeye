#State Space Model Function with autocorrelation
jag.model.SR = function(){  
  
    for (y in 1:nyrs) {
      R[y] ~ dlnorm(mu[y], Tau)
      fit[y] <- log(S[y]) + lnalpha - rk * (beta * S[y])
      e[y] <- log(R[y]) - fit[y]
      }
    mu[1] <- fit[1] + ar1 * phi * e0
    cw[1] ~ dnorm(0, tauw)
    
    for(y in 2:nyrs) {
      cw[y] ~ dnorm(cw[y - 1], tauw)
      mu[y] <- fit[y] + ar1 * phi * e[y - 1]
      }
    lnalpha ~ dunif(min.a, max.a)
    beta ~ dunif(min.b, max.b)
    sigma ~ dunif(0, 2)
    sigmaw ~ dunif(0, 2)
    phi ~ dnorm(0, 5) %_% T(-1, 1)
    e0 ~ dnorm(0, 25)
    Tau <- 1/(sigma * sigma)
    tauw <- 1/(sigmaw * sigmaw)
    lnalpha.c <-lnalpha+0.5*(sigma^2)/(1-phi^2) # mean adjustment
    alpha <- exp(lnalpha)

   }
  