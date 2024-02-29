#' Estimate true causal effect in bivariate Mendelian randomization to adjust for collider bias and weak instrument bias.
#'
#'This function is designed for adjusting collider bias and weak instrument bias to estimate unbiased causal effect between an exposure trait and disease progression trait,
#'using two-step regression and CWLS/MR-RAPS adjustment in bivariate Mendelian randomization.


#' @param data A data frame ofinput dataset. This dataset should contain following columns: \code{xbeta} and \code{xse} represents the effects of exposure trait and its standard errors, \code{dbeta} and \code{dse} represents the effects of disease trait and its standard errors, \code{ybeta} and \code{yse} represents the effects of disease progression trait and its standard errors
#' @param method Adjustment method in each step to correct for collider bias and weak instrument bias. Choose from \code{cwls} and \code{mr.raps}. If MR-RAPS is used, an argument \code{od = TRUE} is needed if invalid instruments are used.
#'
#'
#' @import stats
#'
#'
#' @return An object of class "methodCB" which contains:
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
#' # Find the true causal between exposure and disease progression, using two-step regression with CWLS.
#' TwoSteps(testData, method = "cwls")
#'
#' # Similarly, we can replace the method by MR-RAPS. But note that we are using invalid instrument so over-dispersion needs to be considered.
#' TwoSteps(testData, method = "mr.raps", od = TRUE)
#'
#' @author Siyang Cai, Frank Dudbridge
#'
#' @references
#' Cai S, Hartley A, Mahmoud O, Tilling K, Dudbridge F. Adjusting for collider bias in genetic association studies using instrumental variable methods. Genet Epidemiol. 2022;46:303–16.
#'
#' Zhao, Q, et al. "Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score." The Annals of Statistics 48.3 (2020): 1742-1769.
#'
#' @export
TwoSteps = function(data, method, od = FALSE){
    if(tolower(method) == "cwls"){
        bts = function(data, i = c(1:length(data()))){

            data = data[i,]

            cwls_dy = methodCB(data$dbeta, data$dse, data$ybeta, data$yse, method = "cwls", model = "ivw")

            cwls_reg_beta = methodCB(data$xbeta, data$xse, cwls_dy$ybeta.adj, cwls_dy$yse.adj, method = "cwls", model = "ivw")$b

            return(cwls_reg_beta)
        }

        res = boot::boot(data = data, statistic = bts, R = 100)

        b = res$t0
        b.se = sd(res$t)
    }

    if(tolower(method) == "mr.raps") {
        mr.raps_dy = methodCB(data$dbeta, data$dse, data$ybeta, data$yse, method = "mr.raps", od = od)

        mr.raps_xy = methodCB(data$xbeta, data$xse, mr.raps_dy$ybeta.adj, mr.raps_dy$yse.adj, method = "mr.raps", od = od)

        b = mr.raps_xy$b
        b.se = mr.raps_xy$b.se
    }

    else("Warning: please use CWLS or MR-RAPS to correct for weak instrument bias")

    return(list(b = b, b.se = b.se))
}
