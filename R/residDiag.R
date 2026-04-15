residDiag <- function(mod){
  require(patchwork)
  require(ggplot2)

  if (!(class(mod) %in% c("glm_pois", "glm_pois_zero", "glm_negb", "glm_negb_zero"))){
    stop("Model class not supported")
  }

  # Pearson Residuals Calculations and Plot
  ##########################################
  counts <- mod$model$y
  fit.vals <- as.vector(unlist(mod$fitted.values))
  ris <- ((counts-fit.vals)/sqrt(fit.vals))**2

  p1 <- ggplot(data = data.frame(fit.val=fit.vals,
                                 ri=ris)) +
    geom_smooth(aes(x=fit.val, y=ri),
              color = 'red',
              linetype='dotted') +
    geom_point(aes(x=fit.val, y=ri)) +
    geom_hline(yintercept = 1, linetype = 'dotted') +
    theme_bw() +
    labs(title = expression(Pearson~Residuals~over~mu*"'s"),
         y = expression(r[i]^2),
         x = expression(mu[i]))

  # Randomized Quantized Residuals Calculation and Plot
  ######################################################
  R <- 1000
  rqr <- rep(NA, times=R)
  for (i in 1:R){
    if (class(mod) %in% c("glm_pois", "glm_pois_zero")){
      ai <- ppois(counts[i]-1, lambda = fit.vals[i])
      bi <- ppois(counts[i], lambda = fit.vals[i])
    } else {
      ai <- pnbinom(counts[i]-1, size=mod$theta, mu = fit.vals[i])
      bi <- pnbinom(counts[i], size=mod$theta, mu = fit.vals[i])
    }
    ui <- ai * runif(1) + (bi-ai)
    ui <- max(min(ui, 1-10^(-6)), 10^(-6))
    rqr[i] <- qnorm(ui)
  }

  pearson.ratio <- sum((counts-fit.vals)**2/fit.vals)/mod$df.residual
  p2 <- ggplot(data=data.frame(fit.val=fit.vals,
                               rqrs=rqr)) +
    geom_hline(yintercept=0, linetype="dotted")+
    geom_point(aes(x=fit.val, y=rqrs)) +
    theme_bw() +
    labs(title = "Randomized Quantized Residuals",
          x = bquote(mu),
          y = "RQR")
  p3 <- ggplot(data=data.frame(rqrs=rqr)) +
    stat_qq(aes(sample=rqrs)) +
    stat_qq_line(aes(sample=rqrs)) +
    theme_bw() +
    ggtitle(paste("Dispersion Ratio =", round(pearson.ratio, 4)))

  p1+p2/p3

  diag.message <- ifelse(pearson.ratio>=2, "Overdispersed, try quasi or negative-binomial?",
                     ifelse(pearson.ratio<=0.90, "Underdispersed, try quasi or proceed with caution", "Passed"))
  message(paste("Model fit check:",diag.message))
  message("NOTE: Choose different model structure if wedge shape appears in Pearson Residuals")
}
