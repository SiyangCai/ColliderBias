#' Estimate true causal effect in bivariate Mendelian randomization to adjust for collider bias and weak instrument bias.
#'
#'This function is designed for adjusting collider bias and weak instrument bias to estimate unbiased causal effect between an exposure trait and disease progression trait,
#'using generalised instrumental effect regression and CWBLS adjustment in bivariate Mendelian randomization.


#' @param data A data frame of input dataset. This dataset should contain following columns: \code{xbeta} and \code{xse} represents the effects of exposure trait and its standard errors, \code{dbeta} and \code{dse} represents the effects of disease trait and its standard errors, \code{ybeta} and \code{yse} represents the effects of disease progression trait and its standard errors.
#'
#' @import stats
#'
#'
#' @return An object of class "CWBLS" which contains:
#'   \itemize{
#'     \item{\code{b} The true causal estimate between exposure and disease progression}
#'     \item{\code{b.se} Standard error of \code{b}}
#'   }
#'
#' @examples
#'
#' # Load the test dataset
#' data(testData)
#'
#' # Find the true causal between exposure and disease progression.
#' CWBLS(testData)

#' @author Siyang Cai, Frank Dudbridge
#'
#' @references
#' Cai S, Hartley A, Mahmoud O, Tilling K, Dudbridge F. Adjusting for collider bias in genetic association studies using instrumental variable methods. Genet Epidemiol. 2022;46:303–16.
#'
#'
#' @export
CWBLS = function(data){
  cwbls = function(data, i){
    
    data = data[i,]
  
    w = 1/ data$yse ^ 2
  
    n_i = sum(w * data$dse ^ 2) * sum(w * data$xbeta * data$ybeta)
  
    d_i = sum(w * data$xbeta ^ 2) * sum(w * data$dse ^ 2) + 
      sum(w * data$xse ^ 2) * sum(w * data$dbeta ^ 2) - sum(w * data$xse ^ 2) * sum(w * data$dse ^ 2)
  
    b = (sum(w * data$dbeta ^ 2) * sum(w * data$xbeta * data$ybeta) - sum(w * data$dbeta * data$xbeta) * sum(w * data$dbeta * data$ybeta) - n_i)/
      (sum(w * data$xbeta ^ 2) * sum(w * data$dbeta ^ 2) - (sum(w * data$xbeta * data$dbeta)) ^ 2 - d_i)
  
    return(b)
  }
  
  boot_cwbls = boot::boot(data = data, statistic = cwbls, R = 100, stype = "i")
  
  b = boot_cwbls$t0
  b.se = sd(boot_cwbls$t)
  
  return(list(b = b, b.se = b.se))
}



