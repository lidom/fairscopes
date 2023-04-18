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


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
fair_quantile_EEC_t <- function(alpha, df = NULL, knots, tau, sigma = 1,
                      I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                      alpha_up = alpha*(length(knots)-1), maxIter = 20,
                      tol = alpha / 100){
  # Initialize the u function
  alpha_k = alpha
  ufcns  <- Alg1_gKR_t(alpha = alpha_k, df = df, knots = knots,
                       tau = tau, sigma = sigma, I_weights = I_weights)

  diff <- GeneralKacRice_t(tau = tau,
                           u   = ufcns$u,
                           du  = ufcns$du,
                           df  = df,
                           x   = range(knots),
                           crossings = "up",
                           t0  = 1) - alpha

  niter   = 0
  if(abs(diff) > tol & maxIter != 0){

    alpha_k = alpha_up
    ufcns  <- Alg1_gKR_t(alpha = alpha_k, df = df, knots = knots,
                     tau = tau, sigma = sigma, I_weights = I_weights)

    diff <- GeneralKacRice_t(tau = tau, u = ufcns$u, du = ufcns$du,
                             df = df, x = range(knots),
                             crossings = "up", t0 = 1) - alpha

    if(abs(diff) > tol){
      a = c(alpha, alpha_up)

      while(niter < maxIter & abs(diff) > tol){
        alpha_k = a[1]*0.6 + a[2]*0.4
        ufcns <- Alg1_gKR_t(alpha = alpha_k, df = df, knots = knots,
                        tau = tau, sigma = sigma, I_weights = I_weights)

        diff <- GeneralKacRice_t(tau = tau, u = ufcns$u, du = ufcns$du,
                                 df = df, x = range(knots),
                                 crossings = "up", t0 = 1) - alpha

        if(diff < 0){
          a[1] = alpha_k
        }else{
          a[2] = alpha_k
        }
        niter = niter + 1
      }
    }
  }else{

  }
  return(list(u = ufcns$u, du = ufcns$du, alpha_loc = alpha_k*I_weights, niter = niter))
}
