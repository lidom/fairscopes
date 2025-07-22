#------------------------------------------------------------------------------#
#                                                                              #
#     Bootstrap functions to obtain the quantiles of the maximum of random     #
#     fields                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#      - alg1_KRF()
#      - fair_quantile_KRF()
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' Find an optimal piecewise linear quantile function q to remove
#' conservativeness of standard Kac Rice formula approach for fair
#' thresholds of chi^2 processes.
#'
#' @param alpha threshold to be controlled
#' @param type
#' @param tau  function computing the tau parameter of the process
#' @param I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
#' @param EC.levelset = ">", ...
#' @param ... input for Kac Rice Formula for example for chi^2 the add "df=K"
#' degree of freedom.
#' @return list with elements
#'  \itemize{
#'   \item q Vector containing the fair piecewise linear thresholding function at each x
#'   \item qm Vector containing the offset and the slopes of the fair thresholding function
#'   \item EmpRejections Numeric giving the empirical rejection rate of the fair
#'   thresholding function with respect to the sample.
#' }
#' @export
alg1_KRF <- function(alpha, type, tau, knots,
                        I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                        EC.levelset = ">", df = NULL, tol = alpha/20, ...){
  # If the EC of X>u is computed we need to take lower.tail = FALSE in the point
  # probabilities
  if(EC.levelset == ">"){
    lower.tail = FALSE
    extendInt  = "downX"
    cross      = c("down", "up")
  }else{
    lower.tail = TRUE
    extendInt  = "upX"
    cross      = c("up", "down")
  }
  cross.count = 1

  # Get the KRF for the field of the corresponding type
  #-----------------------------------------------------------------------------
  if(type == "z"){
    KacRice <- function(u, du, x, crossings, t0, EC){
      KacRice_z(tau = tau, u = u, du = du, x = x,
                crossings = crossings, t0 = t0, EC = EC, lower.tail = lower.tail, ...)
    }
  }else if(type == "t"){
    if(is.numeric(df)){
      if(length(df)==1){
        KacRice <- function(u, du, x, crossings, t0, EC){
          KacRice_t(df = df, tau = tau, u = u, du = du, x = x,
                    crossings = crossings, t0 = t0,
                    EC = EC, lower.tail = lower.tail, ...)
        }
      }else{
        stop("The t-field requires setting the degrees of freedom (df) to be a number!")
      }
    }else{
      stop("The t-field requires setting the degrees of freedom (df) to be a number!")
    }
  }else if(type == "elliptic_t"){
    if(is.numeric(df)){
      if(length(df)==1){
        KacRice <- function(u, du, x, crossings, t0, EC){
          KacRice_elliptic_t(df = df, tau = tau, u = u, du = du, x = x,
                             crossings = crossings, t0 = t0,
                             EC = EC, lower.tail = lower.tail, ...)
        }
      }else{
        stop("The elliptic t-field requires setting the degrees of freedom (df) to be a number!")
      }
    }else{
      stop("The elliptic t-field requires setting the degrees of freedom (df) to be a number!")
    }
  }else if(type == "chi2"){
    if(is.numeric(df)){
      if(length(df) == 1){
        KacRice <- function(u, du, x, crossings, t0, EC){
          KacRice_chi2(df = df, tau = tau, u  = u, du = du, x = x,
                       crossings = crossings, t0 = t0, EC = EC, lower.tail = lower.tail)
        }
      }else{
        stop("The chi2-field requires setting the degrees of freedom (df) to be a number!")
      }
    }else{
      stop("The chi2-field requires setting the degrees of freedom (df) to be a number!")
    }
  }else if(type == "F"){
    if(is.numeric(df)){
      if(length(df)==2){
        KacRice <- function(u, du, x, crossings, t0, EC){
          KacRice_F(df = df, tau = tau, u  = u, du = du, x = x,
                    crossings = crossings, t0 = t0, EC = EC, lower.tail = lower.tail, ...)
        }
      }else{
        stop("The F-field requires setting the degrees of freedom (df) to be a numeric 2-vector!")
      }
    }else{
      stop("The F-field requires setting the degrees of freedom (df) to be a numeric 2-vector!")
    }
  }else{
    stop("The type you are requesting is not implmented. Please, consult the help page.")
  }

  # Initialize some paramters
  #-----------------------------------------------------------------------------
  # Get the amount of Intervals
  K     <- length(I_weights)
  # Initialize the parameters for the piecewise linear function on each Interval,
  # first rows are the offsets, second row the slops, columns are the locations/times
  u_fun <- matrix(NA, 2, K)
  # slope = 0 on the first interval
  u_fun[2, 1] <- 0
  # Initialize string for piecewise linear function
  fct_body_u  <- paste0("c_v[1]")
  fct_body_du <- paste0("c_v[1]")

  # Find constant for first interval
  #-----------------------------------------------------------------------------
  find_u0 <- Vectorize(function(u0){
    KacRice(u  = Vectorize(function(y){u0}),
            du = Vectorize(function(y){0}),
            crossings = cross[cross.count],
            x  = c(knots[1], knots[2]), t0 = 2,
            EC = TRUE) - alpha * I_weights[1]
  })
  # get the constant on the first interval
  u_fun[1, 1] <- uniroot(f = find_u0, interval = c(1e-1, 50),
                         extendInt = extendInt, tol = tol)$root

  # Next crossing needs to be of opposite type
  cross.count <- 2

  # Find u on 2:(2*ceiling(K/2)+1) intervals
  #-----------------------------------------------------------------------------
  if(K == 1){
    return(Vectorize(function(t){u_fun[1, 1]}))
  }else if(K != 2){
    for(k in 2:ceiling(K/2)){
      # Even intervals
      #-------------------------------------------------------------------------
      # Initialize start value on the interval
      u_fun[1, 2*k-2] <- u_fun[1, 2*k-3] + u_fun[2, 2*k-3] * (knots[2*k-2] - knots[2*k-3])

      # Function to minimize to get slope
      find_mk <- Vectorize(function(mk){
        KacRice(u  = Vectorize(function(y){u_fun[1, 2*k-2] + mk*(y - knots[2*k-2])}),
                du = Vectorize(function(y){rep(mk, length(y))}),
                crossings = cross[cross.count],
                x  = c(knots[2*k-2], knots[2*k-1]), t0 = 1,
                EC = TRUE) - alpha * I_weights[2*k-2]
      })
      u_fun[2, 2*k-2] <- uniroot(f = find_mk, interval = c(-20, 20),
                                 extendInt = extendInt, tol = tol)$root

      # Next crossing needs to be of opposite type
      cross.count <- 1

      # String for p.linear function
      fct_body_u <- paste0(fct_body_u,
                           "+ c_v[",2*k-2,"]*pmin(pmax(t - knots[",2*k-2,"],0),  knots[",2*k-1,"] - knots[",2*k-2,"])")
      fct_body_du <- paste0(fct_body_du,
                            "+ c_v[",2*k-2,"]*pmax(sign(t - knots[",2*k-2,"]),0)*pmax(-sign(t - knots[",2*k-1,"]),0)")

      # Odd intervals
      #-------------------------------------------------------------------------
      # Initialize first value on the interval
      u_fun[1, 2*k-1] <- u_fun[1, 2*k-2] +  u_fun[2, 2*k-2] * (knots[2*k-1] - knots[2*k-2])

      # Function to minimize to get slope
      find_mk2 <- Vectorize(function(mk){
        KacRice(u  = function(y){u_fun[1, 2*k-1] + mk*(y - knots[2*k-1])},
                du = function(y){rep(mk, length(y))},
                crossings = cross[cross.count],
                x = c(knots[2*k-1], knots[2*k]), t0 = 2,
                EC = TRUE) - alpha * I_weights[2*k-1]
      })
      u_fun[2, 2*k-1] <- uniroot(f = find_mk2, interval = c(-20, 20),
                                 extendInt = extendInt, tol = tol)$root

      # Next crossing needs to be of opposite type
      cross.count <- 2

      # String for p.linear function
      fct_body_u <- paste0(fct_body_u,
                           "+ c_v[",2*k-1,"]*pmin(pmax(t - knots[",2*k-1,"],0),  knots[",2*k,"] - knots[",2*k-1,"])")
      fct_body_du <- paste0(fct_body_du, "+ c_v[",2*k-1,"]*pmax(sign(t - knots[",2*k-1,"]),0)*pmax(-sign(t - knots[",2*k,"]),0)")
    }
  }

  # Find the optimized u for the last interval if there is an even amount of intervals
  #-----------------------------------------------------------------------------
  if(K %% 2 == 0){
    # Initialize value on the interval
    if(K == 2){
      u_fun[1, 2] <- u_fun[1, 1]
    }else{
      u_fun[1, K] <- u_fun[1, K-1] +  u_fun[2, K-1] * (knots[K] - knots[K-1])
    }
    # Function to minimize
    find_mK <- Vectorize(function(mK){
      KacRice(u = function(y){u_fun[1, K] + mK*(y - knots[K])},
              du = function(y){rep(mK, length(y))},
              crossings = cross[cross.count],
              x = c(knots[K], knots[K+1]), t0 = 1,
              EC = TRUE) - alpha * I_weights[K]
    })

    u_fun[2,K] <- uniroot(f = find_mK, interval = c(-10, 10), tol = tol, extendInt = extendInt)$root

    # String for p.linear function
    fct_body_u <- paste0(fct_body_u, "+ c_v[",K,"]*pmin(pmax(t - knots[",K,"],0),  knots[",K+1,"] - knots[",K,"])")
    fct_body_du <- paste0(fct_body_du, "+ c_v[",K,"]*pmax(sign(t - knots[",K,"]),0)*pmax(-sign(t - knots[",K+1,"]),0)")
  }


  # Construct and output the computed function u
  #-----------------------------------------------------------------------------
  coeffs = c(u_fun[1,1], u_fun[2,-1])

  ufun       <- function(t, c_v, knots){}
  body(ufun) <- parse(text = fct_body_u)
  dufun       <- function(t, c_v, knots){}
  body(dufun) <- parse(text = fct_body_du)

  return(list(u = function(t){ ufun(t, c_v = coeffs, knots = knots) },
              du = function(t){ dufun(t, c_v = u_fun[2,], knots = knots) }, c_v =coeffs,
              EC.levelset = EC.levelset))
}

#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams alg1_KRF
#' @param u.type
#' @param alpha_up
#' @param maxIter
#' @param tol
#' @return Standard error under the assumption the data is Gaussian
#' @export
fair_quantile_KRF <- function(alpha, type, tau,
                              knots       = c(0, 1),
                              I_weights   = rep(1 / (length(knots) - 1), length(knots) - 1),
                              u.type      = "linear",
                              alpha_up    = alpha * (length(knots) - 1),
                              EC.levelset = ">",
                              maxIter     = 20,
                              tol         = alpha / 20,
                              df = NULL, ...){
  # If the EC of X>u is computed we need to take lower.tail = FALSE in the point
  # probabilities
  if(EC.levelset == ">"){
    lower.tail = FALSE
    crossings  = "up"
    t0 = 1
  }else{
    lower.tail = TRUE
    crossings  = "down"
    t0 = 1
  }

  # Get the KRF for the field of the corresponding type
  #-----------------------------------------------------------------------------
  type = ifelse(is.null(df) && type %in% c("t", "elliptic_t"), "z", type )

  if(type == "z"){
    KacRice <- function(u, du){
      KacRice_z(tau = tau, u = u, du = du, x = range(knots),
                crossings = crossings, t0 = t0, EC = TRUE, lower.tail = lower.tail)
    }
  }else if(type == "t"){
    if(is.numeric(df)){
      if(length(df) == 1){
        KacRice <- function(u, du){
          KacRice_t(df = df, tau = tau, u = u, du = du, x = range(knots),
                    crossings = crossings, t0 = t0, EC =  TRUE, lower.tail = lower.tail)
        }
      }else{
        stop("The t-field requires setting the degrees of freedom (df) to be a number!")
      }
    }else{
      stop("The t-field requires setting the degrees of freedom (df) to be a number!")
    }
  }else if(type == "elliptic_t"){
    if(is.numeric(df)){
      if(length(df) == 1){
        KacRice <- function(u, du){
          KacRice_elliptic_t(df = df, tau = tau, u = u, du = du, x = range(knots),
                             crossings = crossings, t0 = t0,
                             EC =  TRUE, lower.tail = lower.tail)
        }
      }else{
        stop("The elliptic t-field requires setting the degrees of freedom (df) to be a number!")
      }
    }else{
      stop("The elliptic t-field requires setting the degrees of freedom (df) to be a number!")
    }
  }else if(type == "chi2"){
    if(is.numeric(df)){
      if(length(df) == 1){
        KacRice <- function(u, du){
          KacRice_chi2(df = df, tau = tau, u  = u, du = du, x = range(knots),
                       crossings = crossings, t0 = t0,
                       EC = TRUE, lower.tail = lower.tail)
        }
      }else{
        stop("The chi2-field requires setting the degrees of freedom (df) to be a number!")
      }
    }else{
      stop("The chi2-field requires setting the degrees of freedom (df) to be a number!")
    }
  }else if(type == "F"){
    if(is.numeric(df)){
      if(length(df) == 2){
        KacRice <- function(u, du){
          KacRice_F(df = df, tau = tau, u  = u, du = du, x = range(knots),
                    crossings = crossings, t0 = t0, EC = TRUE, lower.tail = lower.tail)
        }
      }else{
        stop("The F-field requires setting the degrees of freedom (df) to be a numeric 2-vector!")
      }
    }else{
      stop("The F-field requires setting the degrees of freedom (df) to be a numeric 2-vector!")
    }
  }else{
    stop("The type you are requesting is not implmented. Please, consult the help page.")
  }


 # Get the correct algorithm 1 for the class of u
 if(u.type == "linear"){
      alg1 <- function(alpha){
                alg1_KRF(alpha = alpha, type = type, tau = tau,
                         knots = knots, I_weights = I_weights,
                         EC.levelset = EC.levelset, df = df, tol = tol)
      }
 }else{
      alg1 <- function(alpha){
                alg1_KRF_const(alpha = alpha, type = type, tau = tau,
                               knots = knots, I_weights = I_weights,
                               EC.levelset = EC.levelset, df = df, tol = tol)
   }
 }

  # Initialize the u function
  if(maxIter > 0 ){
    diff <- Vectorize(function(alpha_k){
      ufcns  <<- alg1(alpha = alpha_k)
      KacRice(u   = ufcns$u,
              du  = ufcns$du) - alpha
    })
    zz      <- uniroot(f = diff, interval = c(alpha, alpha_up),
                       extendInt = "yes", tol = tol)
    diff    <- zz$f.root
    alpha_k <- zz$root
    niter   <- zz$iter
  }else{
    alpha_k = alpha
    ufcns   <- alg1(alpha = alpha_k)
    diff    <- KacRice(u   = ufcns$u,
                       du  = ufcns$du) - alpha
    niter = 0
  }

  return(list(u = ufcns$u, du = ufcns$du, alpha_loc = alpha_k*I_weights,
              alpha_global = diff + alpha, niter = niter))
}


#' Find an optimal piecewise linear quantile function q to remove
#' conservativeness of standard Kac Rice formula approach for fair
#' thresholds of chi^2 processes.
#'
#' @param alpha threshold to be controlled
#' @param type
#' @param tau  function computing the tau parameter of the process
#' @param I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
#' @param EC.levelset = ">", ...
#' @param ... input for Kac Rice Formula for example for chi^2 the add "df=K"
#' degree of freedom.
#' @return list with elements
#'  \itemize{
#'   \item q Vector containing the fair piecewise linear thresholding function at each x
#'   \item qm Vector containing the offset and the slopes of the fair thresholding function
#'   \item EmpRejections Numeric giving the empirical rejection rate of the fair
#'   thresholding function with respect to the sample.
#' }
#' @export
alg1_KRF <- function(alpha, type, tau, knots,
                     I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                     EC.levelset = ">", basis = NULL, tol = alpha/20, ...){

  # Get derivatives of the basis
  #-----------------------------------------------------------------------------
  ba      <- fd(coef = diag(rep(1, nbasis)), basisobj = basis)
  dbasis  <- deriv.fd(ba, Lfd = 1)
  d2basis <- deriv.fd(ba, Lfd = 2)

  # Get the start value of the minimization problem. We use the constant function
  # which exhibits the correct coverage
  #-----------------------------------------------------------------------------
  # Startwert
  KRFu_up = function(y){
    KacRice_z(
      tau = tau, u = Vectorize(function(x) y),
      du = Vectorize(function(x) 0), x = c(0,1),
      crossings = "up", t0 = 1, sigma = 1,
      lower.tail = FALSE, EC = TRUE) - alpha}
  tt = uniroot(f = KRFu_up, interval = c(0,10))

  x0 <- c(1, rep(tt$root, nbasis))
  q0 <- fd(coef = x0[-1], basisobj = basis)


  # Construct and output the computed function u
  #-----------------------------------------------------------------------------
  res <- nloptr(
    x0 = x0,
    eval_f = eval_f,
    eval_g_eq = eval_g_eq,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA", #  "NLOPT_LD_SLSQP", # "NLOPT_LD_MMA", #
      tol_eq = 1e-4,
      tol_ineq = 1e-4,
      maxeval = 500,
      xtol_rel = 1e-2,
      ftol_rel = 1e-2,
      print_level = 0
    )
  )

  # Construct and output the computed function u
  #-----------------------------------------------------------------------------
  return(list(u = function(t){ ufun(t, c_v = coeffs, knots = knots) },
              du = function(t){ dufun(t, c_v = u_fun[2,], knots = knots) }, c_v =coeffs,
              EC.levelset = EC.levelset))
}
