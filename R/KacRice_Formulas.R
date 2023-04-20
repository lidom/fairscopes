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
              du(t) / sqrt(2*pi*sigma^2) * exp(-u(t)^2 / (2*sigma^2)) *
                pnorm(q = du(t) / sigma / tau(t))
    }

    f3 <-  function(t){
      du(t) / sqrt(2*pi*sigma^2) * exp(-u(t)^2 / (2*sigma^2)) *
        pnorm(q = -du(t) / sigma / tau(t))
    }
  }else{
    nup = df + 1
    afun <- function(t){ sqrt(nu * tau(t)^2 * (1 + u(t)^2 / df) / nup)}

    # Add here the functions for the t-Kac-Rice formula
    f0 <- function(t){
        stats::pt(q = -u(t), df = df)
    }

    f1 <- function(t){
        tau(t) * (1 + u(t)^2 / df + u(t)^2/(df * tau_f(t)^2))^(-df / 2) / (2*pi)
    }

    f2 <- function(t){
      du(t) / (2*pi*tau(t)) * (1 + u(t)^2 / df)^(-df/2 - 1)  *
        gamma(nup/2) * sqrt(nup*pi) * afun(t) / gamma((nup + 1) / 2) *
        stats::pt(q = du(t) / afun(t), df = nup)
    }

    f3 <- function(t){
      du(t) / (2*pi*tau(t)) * (1 + u(t)^2 / df)^(-df/2 - 1)  *
        gamma(nup/2) * sqrt(nup*pi) * afun(t) / gamma((nup + 1) / 2) *
        stats::pt(q = -du(t) / afun(t), df = nup)
    }
  }

  I1 <- f0(-u(x[t0]) / sigma) +
          integrate(f1, lower = x[1], upper = x[length(x)])$value

  if(crossings == "up"){
    return(I1 - integrate(f3, lower = x[1], upper = x[length(x)])$value)
  }else{
    return(I1 + integrate(f2, lower = x[1], upper = x[length(x)])$value)
  }
}
