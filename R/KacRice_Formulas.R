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
    afun <- function(t){ sqrt(df * tau(t)^2 * (1 + u(t)^2 / df) / nup)}

    # Add here the functions for the t-Kac-Rice formula
    f0 <- function(t){
        stats::pt(q = t, df = df)
    }

    f1 <- function(t){
        tau(t) / (2*pi) * (1 + u(t)^2 / df + du(t)^2 / (df * tau(t)^2))^(-df / 2)
    }

    f2 <- function(t){
      du(t) / (2*pi*tau(t)) * (1 + u(t)^2 / df)^(-df/2 - 1)  *
        gamma(nup/2) / gamma((nup + 1) / 2) * sqrt(nup*pi) * afun(t) *
        stats::pt(q = du(t) / afun(t), df = nup)
    }

    f3 <- function(t){
      du(t) / (2*pi*tau(t)) * (1 + u(t)^2 / df)^(-df/2 - 1)  *
        gamma(nup/2) / gamma((nup + 1) / 2) * sqrt(nup*pi) * afun(t) *
        stats::pt(q = -du(t) / afun(t), df = nup)
    }
  }

  If1 <- integrate_save(f1, xlims = c(x[1], x[length(x)]) )

  I1 <- f0(-u(x[t0]) / sigma) + If1

  if(crossings == "up"){
    If3 <- integrate_save(f3, xlims = c(x[1], x[length(x)]) )
    return(I1 - If3 )
  }else{
    If2 <- integrate_save(f2, xlims = c(x[1], x[length(x)]) )
    return(I1 + If2)
  }
}
