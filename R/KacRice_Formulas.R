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

  if(crossings == "down"){
    If2 <- integrate_save(f2, xlims = c(x[1], x[length(x)]) )
    return(I1 + If2 )
  }else{
    If3 <- integrate_save(f3, xlims = c(x[1], x[length(x)]) )
    return(I1 - If3)
  }
}

#' This functions computes the generalization of the Kac-Rice formula
#' from ...
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
KacRice_chisq <- function(tau, u, du, x, df, crossings = "up", t0 = NULL,
                          lower.tail = NULL){
  # Determines at which index of x the marginal probability is computed
  if(is.null(t0)){
    if(crossings == "up"){
      t0 = 1
    }else{
      t0 = length(x)
    }
  }

  # For downcrossings usually the lower excursion is considered for upcrossings
  # usually the upper excursion.
  if(is.null(lower.tail)){
    if(crossings == "up"){
      lower.tail = FALSE
    }else{
      lower.tail = TRUE
    }
  }

  N = df-1

  fac = sqrt(pi)*2^(N/2+1)*gamma((N+1)/2)
  # Integrands in chi^2 Kac Rice formula
  f0 <- Vectorize(function(v) pchisq(q=v, df=df,
                                     lower.tail = lower.tail))

  f1 <- function(t){
    u(t)^((N-2)/2) * exp( -u(t)/2 ) *
          2 * tau(t) * u(t) * exp( -(du(t)/tau(t))^2 / (8*u(t))  ) / fac
  }

  f2 <-  function(t){
    u(t)^((N-2)/2) * exp( -u(t)/2 ) *  sqrt(2*pi)  *
      du(t) * sqrt(u(t)) * pnorm(du(t) / (2 * sqrt(u(t)) * tau(t))) / fac
  }

  f3 <-  function(t){
    u(t)^((N-2)/2) * exp( -u(t)/2 ) * sqrt(2*pi) *
      du(t) * sqrt(u(t)) * pnorm(-du(t) / (2 * sqrt(u(t)) * tau(t))) / fac
  }

  If1 <- integrate_save(f1, xlims = c(x[1], x[length(x)]) )

  if(crossings == "down"){
    If2 <- integrate_save(f2, xlims = c(x[1], x[length(x)]) )
    return( f0(u(x[t0])) + If1 + If2 )
  }else{
    If3 <- integrate_save(f3, xlims = c(x[1], x[length(x)]) )
    return( f0(u(x[t0])) + If1 - If3 )
  }
}


#' This functions computes the generalization of the Kac-Rice formula
#' from ...
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
KacRice_F <- function(tau, u, du, x, df, crossings = "up", t0 = NULL,
                      lower.tail = NULL, add.quantile = FALSE){
  # Determines at which index of x the marginal probability is computed
  if(is.null(t0)){
    if(crossings == "up"){
      t0 = 1
    }else{
      t0 = length(x)
    }
  }

  # For downcrossings usually the lower excursion is considered for upcrossings
  # usually the upper excursion.
  if(is.null(lower.tail)){
    if(crossings == "up"){
      lower.tail = FALSE
    }else{
      lower.tail = TRUE
    }
  }

  N = df[1]
  M = df[2]

  kappa = N/M

  if( N < 200 & M<200){
    fac1 = gamma((N+M-1)/2) / gamma(M/2) / gamma(N/2) / sqrt(pi)
  }else{
    fac1 = (1+M/N)^(N/2-1/2) * (1+N/M)^(M/2-1/2) / sqrt(2*pi) / sqrt(pi)
  }
  if( N < 200 & M<200){
    fac2 = gamma((N+M)/2) / gamma(M/2) / gamma(N/2)
  }else{
    fac2 = (1+M/N)^(N/2-1/2) * (1+N/M)^(M/2-1/2) * sqrt(N+M) / sqrt(2*pi)
  }
  fac1 = fac1
  fac2 = kappa * fac2

  # Integrands in chi^2 Kac Rice formula
  f0 <- Vectorize(function(v) pf( q = v, df1 = N, df2 = M,
                                  lower.tail = lower.tail))

  f1 <- function(t){
    (kappa*u(t))^((N-1)/2) / (1 + kappa*u(t))^((N+M)/2-1) * tau(t) *
      ( 1 + kappa * (du(t)/tau(t))^2 / (4 * u(t) * ( 1 + kappa*u(t))^2) )^(-(N+M-1)/2)
  }

  f2 <-  function(t){
    (kappa*u(t))^((N-2)/2) / ( 1 + kappa*u(t) )^((N+M)/2) * du(t) *
      pt( sqrt(N + M) * sqrt(kappa)* du(t) / ( sqrt(4*u(t)) * (1+kappa*u(t)) * tau(t) ),
          df = N + M, lower.tail = TRUE )
  }

  f3 <-  function(t){
    (kappa*u(t))^((N-2)/2) / ( 1 + kappa*u(t) )^((N+M)/2) * du(t) *
      pt( sqrt(N + M) * sqrt(kappa)* du(t) / ( sqrt(4*u(t)) * (1+kappa*u(t)) * tau(t) ),
          df = N + M, lower.tail = FALSE )
  }

  If1 <- integrate_save(f1, xlims = c(x[1], x[length(x)]) )

  if(add.quantile){
    if(crossings == "down"){
      If2 <- integrate_save(f2, xlims = c(x[1], x[length(x)]) )
      return( f0(u(x[t0])) + fac1*If1 + fac2*If2 )
    }else{
      If3 <- integrate_save(f3, xlims = c(x[1], x[length(x)]) )
      return( f0(u(x[t0])) + fac1*If1 - fac2*If3 )
    }
  }else{
    if(crossings == "down"){
      If2 <- integrate_save(f2, xlims = c(x[1], x[length(x)]) )
      return( fac1*If1 + fac2*If2 )
    }else{
      If3 <- integrate_save(f3, xlims = c(x[1], x[length(x)]) )
      return( fac1*If1 - fac2*If3 )
    }
  }

}
