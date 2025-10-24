#------------------------------------------------------------------------------#
#                                                                              #
#     Bootstrap functions to obtain the quantiles of the maximum of random     #
#     fields                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#      - alg1_KR_chi()
#      - alg1_KR_const()
#      - fair_quantile_EEC_chisq()
#      - alg1_z_DL()
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#

#' Find an optimal piecewise linear quantile function q to remove
#' conservativeness of standard Kac Rice formula approach for fair
#' thresholds of chi^2 processes.
#'
#' @param sample add
#' @return list with elements
#'  \itemize{
#'   \item q Vector containing the fair piecewise linear thresholding function at each x
#'   \item qm Vector containing the offset and the slopes of the fair thresholding function
#'   \item EmpRejections Numeric giving the empirical rejection rate of the fair
#'   thresholding function with respect to the sample.
#' }
#' @export
alg1_KR_chisq <- function(alpha, knots, tau, df,
                          cross.start = "up",
                          I_weights = rep(1/(length(knots) - 1), length(knots) - 1)){
  # Define the even crossing direction
  if(cross.start == "up"){
    cross.even = "down"
  }else if(cross.start == "down"){
    cross.even = "up"
  }else{
    stop("cross.start needs to be either 'up' or 'down'!")
  }

  # Get the amount of Intervals
  K     <- length(I_weights)
  # Initialize the parameters for the piecewise linear function on each Interval
  u_fun <- matrix(NA, 2, K)

  ##############################################################################
  # Find constant for first interval, crossing = can be arbitrary as only matter
  # for non-constant u.
  find_u0 <- Vectorize(function(u0){
    KacRice_chisq(tau = tau,
                     u  = Vectorize(function(y){u0}),
                     du = Vectorize(function(y){0}),
                     x  = c(knots[1], knots[2]),
                     df = df, t0 = 2) - alpha * I_weights[1]
  })
  # slope =0 on the first interval
  u_fun[2, 1] <- 0
  # get the constant on the first interval
  u_fun[1, 1] <- uniroot(f = find_u0, interval = c(1e-1, 50),
                         extendInt = "downX")$root

  # Initialize string for piecewise linear function
  fct_body_u  <- paste0("c_v[1]")
  fct_body_du <- paste0("c_v[1]")

  ##############################################################################
  # Find constant for first interval
  if(K == 1){
    return(Vectorize(function(t){u_fun[1, 1]}))
  }else if(K != 2){
    for(k in 2:ceiling(K/2)){
      #---- even intervals
      # Initialize start value on the interal
      u_fun[1, 2*k-2] <- u_fun[1, 2*k-3] + u_fun[2, 2*k-3] * (knots[2*k-2] - knots[2*k-3])

      # Function to minimize to get slope
      find_mk <- Vectorize(function(mk){
        KacRice_chisq(tau,
                      u  = Vectorize(function(y){u_fun[1, 2*k-2] + mk*(y - knots[2*k-2])}),
                      du = Vectorize(function(y){mk}),
                      x  = c(knots[2*k-2], knots[2*k-1]),
                      df = df, crossings = cross.even, t0 = 1) - alpha * I_weights[2*k-2]
      })
      u_fun[2, 2*k-2] <- uniroot(f = find_mk, interval = c(-20, 20), extendInt = "downX")$root

      # String for p.linear function
      fct_body_u <- paste0(fct_body_u,
                           "+ c_v[",2*k-2,"]*pmin(pmax(t - knots[",2*k-2,"],0),  knots[",2*k-1,"] - knots[",2*k-2,"])")
      fct_body_du <- paste0(fct_body_du,
                            "+ c_v[",2*k-2,"]*pmax(sign(t - knots[",2*k-2,"]),0)*pmax(-sign(t - knots[",2*k-1,"]),0)")

      #---- odd intervals
      # Initialize first value on the interval
      u_fun[1, 2*k-1] <- u_fun[1, 2*k-2] +  u_fun[2, 2*k-2] * (knots[2*k-1] - knots[2*k-2])

      # Function to minimize to get slope
      find_mk2 <- Vectorize(function(mk){
        KacRice_chisq(tau,
                         u  = function(y){u_fun[1, 2*k-1] + mk*(y - knots[2*k-1])},
                         du = function(y){rep(mk, length(y))}, x = c(knots[2*k-1], knots[2*k]),
                         df = df, crossings = "down", t0=2) - alpha * I_weights[2*k-1]
      })
      u_fun[2, 2*k-1] <- uniroot(f = find_mk2, interval = c(-20, 20), extendInt = "downX")$root

      # String for p.linear function
      fct_body_u <- paste0(fct_body_u,
                           "+ c_v[",2*k-1,"]*pmin(pmax(t - knots[",2*k-1,"],0),  knots[",2*k,"] - knots[",2*k-1,"])")
      fct_body_du <- paste0(fct_body_du, "+ c_v[",2*k-1,"]*pmax(sign(t - knots[",2*k-1,"]),0)*pmax(-sign(t - knots[",2*k,"]),0)")
    }
  }

  # Find the optimized u for the last interval if there is an even amount of intervals
  if(K %% 2 == 0){
    # Initialize value on the interval
    if(K == 2){
      u_fun[1, 2] <- u_fun[1, 1]
    }else{
      u_fun[1, K] <- u_fun[1, K-1] +  u_fun[2, K-1] * (knots[K] - knots[K-1])
    }
    # Function to minimize
    find_mK <- Vectorize(function(mK){
      KacRice_chisq(tau, u = function(y){u_fun[1, K] + mK*(y - knots[K])},
                       du = function(y){rep(mK, length(y))}, x = c(knots[K], knots[K+1]),
                       df = df, crossings = "up") - alpha * I_weights[K]
    })
    u_fun[2,K] <- uniroot(f = find_mK, interval = c(-20, 20))$root

    # String for p.linear function
    fct_body_u <- paste0(fct_body_u, "+ c_v[",K,"]*pmin(pmax(t - knots[",K,"],0),  knots[",K+1,"] - knots[",K,"])")
    fct_body_du <- paste0(fct_body_du, "+ c_v[",K,"]*pmax(sign(t - knots[",K,"]),0)*pmax(-sign(t - knots[",K+1,"]),0)")
  }

  coeffs = c(u_fun[1,1], u_fun[2,-1])

  ufun       <- function(t, c_v, knots){}
  body(ufun) <- parse(text = fct_body_u)
  dufun       <- function(t, c_v, knots){}
  body(dufun) <- parse(text = fct_body_du)

  return(list(u = function(t){ ufun(t, c_v = coeffs, knots = knots) },
              du = function(t){ dufun(t, c_v = u_fun[2,], knots = knots) }, c_v =coeffs))
}


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @param cross.start can be "up" or "down". Do we count up or downcrossings on
#' the first interval.
#' @return Standard error under the assumption the data is Gaussian
#' @export
fair_quantile_EEC_chisq <- function(alpha, tau, x = seq(0, 1, length.out = 2), df,
                                crit.set    = list(minus = rep(T, length(x)),
                                                   plus  = rep(T, length(x))),
                                knots       = range(x),
                                I_weights   = rep(1/(length(knots) - 1), length(knots) - 1),
                                type        = "linear",
                                alpha_up    = alpha*(length(knots)-1),
                                crossings   = "up",
                                maxIter     = 20,
                                tol         = alpha / 20){
  # Get the correct algorithm 1 for the class of u and whether SCoPES or SCBs
  # are computed
  if(all(crit.set$minus) | all(crit.set$plus)){
    if(type == "linear"){
      alg1 <- alg1_KR_chisq
    }else{
      alg1 <- alg1_KR_chisq_const
    }
  }else{
    # Not coded
    #---------------------------------------------------------------------------
    knots     <- adapt_knots(knots, crit.set)
    I_weights <- adapt_Iweights(I_weights, crit.set)

    alg1 <- alg1_KR_t_pw
    #---------------------------------------------------------------------------
  }


  # Initialize the u function
  alpha_k = alpha

  ufcns  <- alg1(alpha = alpha_k, tau = tau, df = df,
                 knots = knots, crossings = crossings, I_weights = I_weights)

  diff <- KacRice_chisq(tau = tau,
                        u   = ufcns$u,
                        du  = ufcns$du,
                        df  = df,
                        x   = range(knots),
                        crossings = crossings,
                        t0  = 1) - alpha

  niter   = 0
  if(abs(diff) > tol & maxIter != 0){

    alpha_k = alpha_up
    ufcns  <- alg1(alpha = alpha_k, tau = tau, df = df, crossings = crossings,
                   knots = knots, I_weights = I_weights)

    diff <- KacRice_chisq(tau = tau, u = ufcns$u, du = ufcns$du,
                          df  = df,  x = range(knots),
                          crossings = crossings, t0 = 1) - alpha

    if(abs(diff) > tol){
      a = c(alpha, alpha_up)

      while(niter < maxIter & abs(diff) > tol){
        if(niter < 3){
          alpha_k = a[1]*0.8 + a[2]*0.2
        }else{
          alpha_k = mean(a)
        }

        ufcns <- alg1(alpha = alpha_k, tau = tau, df = df, crossings = crossings,
                      knots = knots, I_weights = I_weights)

        diff <- KacRice_chisq(tau = tau, u = ufcns$u, du = ufcns$du,
                              df = df, x = range(knots),
                              crossings = crossings, t0 = 1) - alpha

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
  return(list(u = ufcns$u, du = ufcns$du, alpha_loc = alpha_k*I_weights,
              alpha_global = diff+alpha, niter = niter))
}
