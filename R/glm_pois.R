glm_pois <- function(data,
                     formula,
                     offset = log(1)){

  #Parameter initliazations
  ##########################
  par <- model.frame(formula, data=data)
  log.lambdas <- model.response(par)
  X <- model.matrix(formula, data=par)
  betas <- matrix(0, nrow=ncol(X),ncol=1)

  #IWLS algorithm model fit
  ##########################
  repeat{
    eta <- offset + X %*% betas
    pred.means <- exp(eta)
    W <- diag(pred.means)
    z <- eta+(log.lambdas-pred.means)/pred.means
    inv.Hess <- solve(t(X)%*%W%*%X)
    betas.new <- inv.Hess%*%t(X)%*%W%*%z
    ss <- sum((betas.new-betas)**2)
    betas <- betas.new
    std.error <- sqrt(diag(inv.Hess))
    if(ss < 1e-6){
      break
    }
  }

  fit.dat <- as.data.frame(cbind(betas, std.error))

  # Coefficient Pvals and Confidence Intervals
  #############################################
  p.vals <- rep(NA, times = nrow(fit.dat))
  asymp.CI.lower <- rep(NA, times = nrow(fit.dat))
  asymp.CI.higher <- rep(NA, times = nrow(fit.dat))
  SEs <- sqrt(inv.Hess[i,i])
  for(i in 1:nrow(fit.dat)){
    SE <- SEs[i,i]
    Z <- fit.dat$betas[i]/SE
    p.vals[i] <- 2*(1-pnorm(abs(Z)))
    crit <- qnorm(0.975)
    asymp.CI.lower[i] <- fit.dat$betas[i]-crit*SE
    asymp.CI.higher[i] <- fit.dat$betas[i]+crit*SE
  }

  fit.dat <- as.data.frame(cbind(fit.dat, asymp.CI.lower, asymp.CI.higher, p.vals))

  return(fit.dat)
}
