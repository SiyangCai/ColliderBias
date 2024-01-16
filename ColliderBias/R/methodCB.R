#' Adjust association statistics for collider bias and weak instrument bias
#'
#' Given effect sizes and standard errors for predictors of an index trait and a subsequent trait,
#' this function adjusts the statistics for the subsequent trait for selection bias through the index trait.
#'
#' Effect sizes are on a linear scale, so could be the coefficients from linear regression, or log odds ratios, or log hazard ratios.
#' Effects on the subsequent trait are regressed on the effects on the index trait.
#' By default, the regression is weighted by the inverse variances of the subsequent trait effects.
#' The regression is adjusted for sampling variation in the index trait effects,
#' and the residuals then used to obtain adjusted effect sizes and standard errors for the subsequent trait.
#'
#' The regression should be performed on a subset of predictors that are independent.
#' In the context of a genome-wide association study, these would be LD-pruned SNPs.
#' In terms of the input parameters, the regression command is \code{lm(ybeta[prune]~xbeta[prune],weights=1/yse[prune]^2)}.
#'
#' The effects in \code{xbeta} and \code{ybeta} should be aligned for the same variables
#' and the same direction prior to running \code{indexevent}.
#'
#' The default value of \code{B} is 100 to get a quick result, but higher values are recommended, eg 10000.

#' @param xbeta Vector of effects on the index trait.
#' @param xse Vector of standard errors of \code{xbeta}.
#' @param ybeta Vector of effects on the subsequent trait.
#' @param yse Vector of standard errors of \code{ybeta}.
#' @param prune Vector containing the indices of an approximately independent subset of the predictors in \code{xbeta} and \code{ybeta}. If unspecified, all predictors will be used.
#' @param method Method to adjust for regression dilution in the regression of \code{ybeta[prune]} on \code{xbeta[prune]}.
#' "CWLS" applies corrected weighted least squares estimator, "Fast" applies the fast version of CWLS, "SH" applies Slope Hunter and "mr.raps" applies MR-RAPS.
#' @param weighted Weights to be used in raw regression, Fast SIMEX and Hedges-Olkin adjustment. Default is using "1" as the first order weight, while "0" indicates unweighted model and "2" applies second order weight.
#' @param model Rgression model to be used in raw regression, Fast SIMEX and Hedges-Olkin adjustment. Default is "IVW" and could be replaced by "Egger" to use MR-Egger model.
#' @param B Number of simulations performed in each stage of the Fast SIMEX adjustment.
#' @param lambda Vector of lambdas for which the Fast SIMEX simulations are performed.
#' @param seed Random number seed for the Fast SIMEX adjustment.
#' @param od Whether to use "over.dispersion" in MR-RAPS.
#'
#'
#' @import stats
#'
#' @return An object of class "methodCB" which contains:
#'   \itemize{
#'     \item{\code{ybeta.adj} Adjusted effects on the subsequent trait}
#'     \item{\code{yse.adj} Adjusted standard errors of \code{ybeta.adj}}
#'     \item{\code{ychisq.adj} Chi-square statistics for \code{(ybeta.adj/yse.adj)^2}}
#'     \item{\code{yp.adj} P-values for \code{ychisq.adj} on 1df}
#'     \item{\code{b} Coefficient of the regression of \code{ybeta[prune]} on \code{xbeta[prune]}, after correction for regression dilution}
#'     \item{\code{b.se} Standard error of \code{b}}
#'     \item{\code{b.ci} Lower and upper confidence limits for \code{b}}
#'     \item{\code{b.raw} Regression coefficient without correction for regression dilution}
#'     \item{\code{b.raw.se} Standard error of \code{b.raw}}
#'   }
#'
#' @author Siyang Cai, Frank Dudbridge
#'
#' @references
#' Cai S, Hartley A, Mahmoud O, Tilling K, Dudbridge F. Adjusting for collider bias in genetic association studies using instrumental variable methods. Genet Epidemiol. 2022;46:303–16.
#'
#' Mahmoud O, Dudbridge F, Smith G D, et al. Slope-Hunter: A robust method for index-event bias correction in genome-wide association studies of subsequent traits[J]. bioRxiv, 2020.
#'
#' Zhao, Q, et al. "Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score." The Annals of Statistics 48.3 (2020): 1742-1769.
#'
#' @export
methodCB = function (xbeta,
                      xse,
                      ybeta,
                      yse,
                      prune = NULL,
                      method,
                      weighted = "1",
                      model = "ivw",
                      B = 100,
                      lambda = seq(0.25, 5, 0.25),
                      od = FALSE,
                      seed = 0)

{
  if (is.null(prune)){
    prune = 1:length(xbeta)}

  xbetaprune = xbeta[prune]
  xseprune = xse[prune]
  ybetaprune = ybeta[prune]
  yseprune = yse[prune]

  if (weighted == "1"){
    weight = 1 / yseprune ^ 2
  }
  if (weighted == "2"){
    weight = 1 / (yseprune ^ 2 + (ybetaprune / xbetaprune) ^ 2 * xseprune ^ 2)
  }
  if (weighted == "0"){
    weight = rep(1, length(prune))}

  if (tolower(model) == "egger"){
    fit = lm(ybetaprune ~ xbetaprune, weights = weight)
    fits = summary(lm(ybetaprune ~ xbetaprune, weights = weight))
    b.raw = as.numeric(fits$coef[2])
    b.raw.se = as.numeric(fits$coef[4])
  }

  if (tolower(model) == "ivw"){
    fit = lm(ybetaprune ~ xbetaprune - 1, weights = weight)
    fits = summary(lm(ybetaprune ~ xbetaprune - 1, weights = weight))
    b.raw = as.numeric(fits$coef[1])
    b.raw.se = as.numeric(fits$coef[2])
  }

  if (tolower(method) == "cwls") {
    if(tolower(model) == "egger"){
      ho = (sum(weight, na.rm = TRUE) * sum(weight * xbetaprune ^ 2, na.rm = TRUE) - (sum(weight * xbetaprune, na.rm = TRUE)) ^ 2) /
        (sum(weight, na.rm = TRUE) * sum(weight * xbetaprune ^ 2, na.rm = TRUE) - (sum(weight * xbetaprune, na.rm = TRUE)) ^ 2 - sum(weight) * sum(weight * xseprune ^ 2, na.rm = TRUE))
    }

    if (tolower(model) == "ivw"){
      ho = sum(weight * xbetaprune ^ 2, na.rm = TRUE) /
        (sum(weight * xbetaprune ^ 2, na.rm = TRUE) - sum(weight * xseprune ^ 2, na.rm = TRUE))
    }

    b = b.raw * ho
    slope.ho = b
    b.se = b.raw.se * ho
    b.ci = c(b - qnorm(0.975) * b.se, b + qnorm(0.975) * b.se)
    simex.estimates = NULL
  }


  if(tolower(method) == "fast"){
    if(tolower(model) == "egger"){
      simex.estimates = fit$coef[2]
      simex.variance.sandwich = sandwich::vcovHC(fit)[2,2]

      sum.weight.weight <- sum(weight) * weight
      sum.ybeta <- sum(weight * ybetaprune)

      ho = (sum(weight, na.rm = TRUE) * sum(weight * xbetaprune ^ 2, na.rm = TRUE) - (sum(weight * xbetaprune, na.rm = TRUE)) ^ 2) /
        (sum(weight, na.rm = TRUE) * sum(weight * xbetaprune ^ 2, na.rm = TRUE) - (sum(weight * xbetaprune, na.rm = TRUE)) ^ 2 - sum(weight) * sum(weight * xseprune ^ 2, na.rm = TRUE))
      slope.ho = b.raw * ho
      m = length(xbetaprune)

      progress = txtProgressBar(max = length(lambda), width = 10, style = 3)

      set.seed(seed)

      for(l in 1:length(lambda)) {
        # simulate matrix of xbeta + sampling errors
        simex.errors <- matrix(rnorm(m*B, mean= xbetaprune, sd = xseprune * sqrt(lambda[l])), nrow = m, ncol = B)
        weight.simex.errors <- weight %*% simex.errors

        simex.numer <- (sum.weight.weight * ybetaprune) %*% simex.errors - weight.simex.errors * sum.ybeta
        simex.denom <- sum.weight.weight %*% simex.errors^2 - weight.simex.errors^2
        simexcoef <- simex.numer / simex.denom

        # take their mean
        simex.estimates = c(simex.estimates, mean(simexcoef))

        setTxtProgressBar(progress, l)
      }

      set.seed(seed)

      for(l in 1:length(lambda)){
        simexdata = rnorm(length(xbetaprune), mean = xbetaprune, sd = xseprune * sqrt(lambda[l]))
        simexfit = lm(ybetaprune ~ simexdata, weights = weight)
        svar = sandwich::vcovHC(simexfit)[2,2]

        simex.variance.sandwich = c(simex.variance.sandwich, svar/B)
        setTxtProgressBar(progress, l)
      }


      simex.estimates <- cbind(c(0,lambda), simex.estimates, simex.variance.sandwich)
      colnames(simex.estimates) <- c("Lambda", "Coefficient", "Variance")
      rownames(simex.estimates) <- NULL

      simexMLE.fast = nlm(simexllhdbivariate, c(slope.ho ,log(ho-1)), simex.estimates = simex.estimates)


      b = simexMLE.fast$estimate[1]

      h = numDeriv::hessian(simexllhdbivariate, simexMLE.fast$estimate, simex.estimates = simex.estimates)
      h.inv = solve(h)
      b.se = sqrt(h.inv[1,1])
      b.ci = c(b - qnorm(0.975) * b.se, b + qnorm(0.975) * b.se)
    }

    if(tolower(model) == "ivw"){
      simex.estimates = fit$coef
      simex.variance.sandwich = sandwich::vcovHC(fit)
      m = length(xbetaprune)

      ho = sum(weight * xbetaprune ^ 2, na.rm = TRUE) /
        (sum(weight * xbetaprune ^ 2, na.rm = TRUE) - sum(weight * xseprune ^ 2, na.rm = TRUE))
      slope.ho = b.raw * ho

      progress = txtProgressBar(max = length(lambda), width = 10, style = 3)

      set.seed(seed)
      for(l in 1:length(lambda)) {
        # simulate matrix of xbeta + sampling errors
        simex.errors <- matrix(rnorm(m*B, mean= xbetaprune, sd = xseprune * sqrt(lambda[l])), nrow = m, ncol = B)

        simex.numer <- (weight * ybetaprune) %*% simex.errors
        simex.denom <- weight %*% simex.errors ^ 2
        simexcoef <- simex.numer / simex.denom

        # take their mean
        simex.estimates = c(simex.estimates, mean(simexcoef))

        setTxtProgressBar(progress, l)
      }

      set.seed(seed)
      for(l in 1:length(lambda)){
        simexdata = rnorm(length(xbetaprune), mean = xbetaprune, sd = xseprune * sqrt(lambda[l]))
        simexfit = lm(ybetaprune ~ simexdata - 1, weights = weight)
        svar = sandwich::vcovHC(simexfit)[1,1]

        simex.variance.sandwich = c(simex.variance.sandwich, svar/B)

        setTxtProgressBar(progress, l)
      }


      simex.estimates <- cbind(c(0,lambda), simex.estimates, simex.variance.sandwich)
      colnames(simex.estimates) <- c("Lambda", "Coefficient", "Variance")
      rownames(simex.estimates) <- NULL

      simexMLE.fast <- nlm(simexllhdbivariate, c(slope.ho, log(ho - 1)), simex.estimates = simex.estimates)
      b = simexMLE.fast$estimate[1] #-0.3354914

      h <- numDeriv::hessian(simexllhdbivariate, simexMLE.fast$estimate, simex.estimates = simex.estimates)
      h.inv <- solve(h)
      b.se = sqrt(h.inv[1,1])
      b.ci = c(b - qnorm(0.975) * b.se, b + qnorm(0.975) * b.se)
    }
  }

  if(tolower(method) == "sh"){
    SHData = data.frame(xbetaprune, xseprune, ybetaprune, yseprune)
    names(SHData) = c("xbeta","xse", "ybeta", "yse")
    SHresults = SlopeHunter::hunt(dat = SHData,
                                  xbeta_col="xbeta", xse_col="xse",
                                  ybeta_col="ybeta", yse_col="yse",
                                  xp_thresh = 0.05, Bootstrapping = TRUE, show_adjustments = TRUE, Plot = FALSE, seed=2019)
    b = SHresults$b
    b.se = SHresults$bse
    b.ci = SHresults$b_CI

  }

  if(tolower(method) == "mr.raps"){
    MRRAPS.results = mr.raps::mr.raps(xbetaprune, ybetaprune, xseprune, yseprune, over.dispersion = od)
    b = MRRAPS.results$beta.hat
    b.se = MRRAPS.results$beta.se
    b.ci = c(b - qnorm(0.975) * b.se, b + qnorm(0.975) * b.se)

  }


  ybeta.adj = rep(0, length(ybeta))
  ybeta.adj = ybeta - b * xbeta
  yse.adj = sqrt(yse ^ 2 + b ^ 2 * xse ^ 2)
  ychisq.adj = (ybeta.adj / yse.adj) ^ 2
  yp.adj = pchisq(ychisq.adj, 1, lower = F)
  results = list(
    ybeta.adj = ybeta.adj,
    yse.adj = yse.adj,
    ychisq.adj = ychisq.adj,
    yp.adj = yp.adj,
    b = b,
    b.se = b.se,
    b.ci = b.ci,
    b.raw = b.raw,
    b.raw.se = b.raw.se
  )
  class(results) = ("MRmethods")
  results
}
