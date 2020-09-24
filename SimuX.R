set.seed(1)
library(MASS)
###Simulation data
SimuX <- function(Sample_size, Number_predictors, rho){
  n <- Sample_size #sample size
  d <- Number_predictors #number of predictors
  r <- rho
  M <- diag(1,d)
  
  for (i in 1:d)
  {
    for (j in 1:i)
    {
      M[j,i] <- r^{i-j}
      M[i,j] <- M[j,i]
    }
  }
  return(mvrnorm(n,rep(0,d),M)) 
}

