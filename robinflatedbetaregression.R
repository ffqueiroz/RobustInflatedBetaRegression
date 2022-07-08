

fit.betainflated <- function(y, S, X, Z,  alphaoptimal = TRUE, alpha0 = 0, weights = FALSE,
                             linkmu="logit", linkphi="log", linkvartheta = "logit"){
  source("auxiliaryfunctions_RobInfBetaReg.R")
  index0 <- y!=0 # y \in (0, 1)
  
  y.1 <- y[index0]
  X.1 <- X[index0, ]
  Z.1 <- as.matrix(Z[index0, ])
  
  y.2 <- ifelse(index0, 0, 1) # I(y=0)
  
  fit.LSMLE <- LSMLE_BETA(y.1, X.1, Z.1, alphaoptimal = alphaoptimal, alpha0 = alpha0, weights = weights,
                          linkmu=linkmu, linkphi=linkphi)
  fit.LMDPDE <- LMDPDE_BETA(y.1, X.1, Z.1, alphaoptimal = alphaoptimal, alpha0 = alpha0, weights = weights,
                            linkmu=linkmu, linkphi=linkphi)
  fit.SMLE  <- SMLE_BETA(y.1, X.1, Z.1, qoptimal = alphaoptimal, q0 = 1-alpha0, weights = weights,
                         linkmu=linkmu, linkphi=linkphi)
  fit.MDPDE <- MDPDE_BETA(y.1, X.1, Z.1, qoptimal = alphaoptimal, q0 = 1-alpha0, weights = weights,
                          linkmu=linkmu, linkphi=linkphi)
  fit.MDPDE.disc <- MDPDE_BIN(y.2, S, alphaoptimal = alphaoptimal, alpha0 = alpha0, weights = weights,
                              linkvartheta = linkvartheta)
  
  Upsilon.M_SE  <- c(fit.MDPDE.disc$kappa, fit.SMLE$beta, fit.SMLE$gama)
  Upsilon.M_LSE <- c(fit.MDPDE.disc$kappa, fit.LSMLE$beta, fit.LSMLE$gama)
  Upsilon.M_ME  <- c(fit.MDPDE.disc$kappa, fit.MDPDE$beta, fit.MDPDE$gama)
  Upsilon.M_LME <- c(fit.MDPDE.disc$kappa, fit.LMDPDE$beta, fit.LMDPDE$gama)
  
  alpha.LSMLE  <- fit.LSMLE$alpha_value
  alpha.SMLE   <- 1 - fit.SMLE$q_value
  alpha.MDPDE  <- 1 - fit.MDPDE$q_value
  alpha.LMDPDE <- fit.LMDPDE$alpha_value
  
  se.M_SE  <- c(fit.MDPDE.disc$se_kappa, vcov_RobInfBeta(Upsilon.M_SE, "M_SE", alpha.SMLE, S = S, X = X, Z = Z, y = y))
  se.M_LSE <- c(fit.MDPDE.disc$se_kappa, vcov_RobInfBeta(Upsilon.M_LSE, "M_LSE", alpha.LSMLE, S = S, X = X, Z = Z, y = y))
  se.M_ME  <- c(fit.MDPDE.disc$se_kappa, vcov_RobInfBeta(Upsilon.M_ME, "M_ME", alpha.MDPDE, S = S, X = X, Z = Z, y = y))
  se.M_LME <- c(fit.MDPDE.disc$se_kappa, vcov_RobInfBeta(Upsilon.M_LME, "M_LME", alpha.LMDPDE, S = S, X = X, Z = Z, y = y))
  
  estimates <- cbind(Upsilon.M_SE, se.M_SE, 
                     Upsilon.M_LSE, se.M_LSE,
                     Upsilon.M_ME, se.M_ME,
                     Upsilon.M_LME, se.M_LME)  
  
  cat("\n\n ######## Results M-SE ######## \n")
  cat("Estimates\n")
  print(structure(round(as.vector(Upsilon.M_SE), digits = 4)))
  cat("Standard errors\n")
  print(structure(round(as.vector(se.M_SE), digits = 4)))
  cat("Optimum alpha MDPDE (discrete part) = ", fit.MDPDE.disc$alpha_value, "\n")
  cat("Optimum alpha SMLE (continuous part) = ", alpha.SMLE, "\n")
  
  cat("\n\n ######## Results M-ME ######## \n")
  cat("Estimates\n")
  print(structure(round(as.vector(Upsilon.M_ME), digits = 4)))
  cat("Standard errors\n")
  print(structure(round(as.vector(se.M_ME), digits = 4)))
  cat("Optimum alpha MDPDE (discrete part) = ", fit.MDPDE.disc$alpha_value, "\n")
  cat("Optimum alpha MDPDE (continuous part) = ", alpha.MDPDE, "\n")
  
  return(list(estimates = estimates,
              MDPDE.disc = fit.MDPDE.disc,
              LSMLE.cont = fit.LSMLE,
              LMDPDE.cont = fit.LMDPDE,
              SMLE.cont = fit.SMLE,
              MDPDE.cont = fit.MDPDE))
}
