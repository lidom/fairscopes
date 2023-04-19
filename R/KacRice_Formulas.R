#------------------------------------------------------------------------------#
#                                                                              #
#     Generalized Kac-Rice Formulas and Fair EECs
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
# - GeneralKacRice_t()
# - fair_quantile_EEC_t()
#------------------------------------------------------------------------------#
# Developer notes:
# - Fix the documentation and man pages
#
#------------------------------------------------------------------------------#
#' This functions computes the generalization of the Kac-Rice formula
#' from ...
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
GeneralKacRice_t <- function(tau, u, du, x, df = NULL,
                  sigma = 1, crossings = "up", t0 = NULL){
  # Determines at which index of x the marginal probability is computed
  if(is.null(t0)){
    if(crossings == "up"){
      t0 = 1
    }else{
      t0 = length(x)
    }
  }

  # Integrands in Dominik's Kac Rice formula
  if(is.null(df)){
    f0 <- pnorm

    f1 <- function(t){
      tau(t) / (2*pi) * exp( -(u(t)^2 + (du(t)/tau(t))^2) / (2*sigma^2) )
    }

    f2 <-  function(t){
      du(t) / sqrt(2*pi*sigma^2) *
        exp(-u(t)^2 / (2*sigma^2)) *
        pnorm(q = du(t) / sigma / tau(t))
    }

    f3 <-  function(t){
      du(t) / sqrt(2*pi*sigma^2) *
        exp(-u(t)^2 / (2*sigma^2)) *
        pnorm(q = -du(t) / sigma / tau(t))
    }
  }else{
    # Add here the functions for the t-Kac-Rice formula
  }

  I1 <- f0(-u(x[t0]) / sigma) +
          integrate(f1, lower = x[1], upper = x[length(x)])$value

  if(crossings == "up"){
    return(I1 - integrate(f3, lower = x[1], upper = x[length(x)])$value)
  }else{
    return(I1 + integrate(f2, lower = x[1], upper = x[length(x)])$value)
  }
}
