#------------------------------------------------------------------------------#
#                                                                              #
#     Bootstrap functions to obtain the quantiles of the maximum of random     #
#     fields                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#      - fair_Bootstrap()
#      - fair_quantile_boot()
#      - quantile_KRF()
#      - fair_quantile_KRF()
#      - alg1_KRF_t()
#      - alg1_KRF_const()
#      - alg1_z_DL()
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------
# Bootstrap Algorithms
#-------------------------------------------------------------------------------
#' Estimates from a sample of random functions (for example a bootstrap sample)
#' a Fair thresholding function q
#'
#' @param sample array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param x vector of length K of locations at which the sample is observed.
#' @param alpha numeric the targeted upper quantile of the maxiumum of the
#'   absolute value of the random field. Default is 0.95.
#' @param fair a vector partitioning the vector x into regions on which the
#'             on which the rejection should be fair. First element must be
#'             x[1] and last x[length(x)]
#' @return list with elements
#'  \itemize{
#'   \item q Vector containing the fair piecewise linear thresholding function at each x
#'   \item qm Vector containing the offset and the slopes of the fair thresholding function
#'   \item EmpRejections Numeric giving the empirical rejection rate of the fair
#'   thresholding function with respect to the sample.
#' }
#' @export
fair_Bootstrap <- function(alpha, samples, x,
                           crit.set  = list(minus = rep(T, length(x)),
                                            plus  = rep(T, length(x))),
                           knots     = range(x),
                           I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                           type      = "linear",
                           subI      = NULL,
                           inter     = NULL ){
  #
  if(is.null(subI) || is.null(inter)){
    s     = sub.intervals(x, knots, crit.set)
    subI  = s$subI
    inter = s$inter
  }

  # Simplify notation
  Sminus = crit.set$minus
  Splus  = crit.set$plus

  # Spread the parts where no probability is expected fairly
  # to the other partitions
  if(!all(inter)){
    I_weights = I_weights + sum(I_weights[!inter]) / sum(inter)
    I_weights[!inter] = 0
  }

  dimS = dim(samples)
  #  if(fair.type = "linear"){
  samples_minus = -samples
  samples_minus[!Sminus,] <- -Inf
  samples_plus = samples
  samples_plus[!Splus,]   <- -Inf
  #}

  subI_prob <- function(k, q){
    if(is.null(Sminus)){
      low.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_minus     = intersect(which(Sminus), subI[[k]])
      if(length(ind_minus) == 0){
        low.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_minus) == 1){
        low.excursion = -samples[ind_minus,] - q[ind_minus] >= 0
      }else{
        low.excursion = apply(-samples[ind_minus,] - q[ind_minus], 2, max) >= 0
      }
    }

    if(is.null(Splus)){
      up.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_plus     = intersect(which(Splus), subI[[k]])
      if(length(ind_plus) == 0){
        up.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_plus) == 1){
        up.excursion = samples[ind_plus,] - q[ind_plus] >= 0
      }else{
        up.excursion = apply(samples[ind_plus,] - q[ind_plus], 2, max) >= 0
      }
    }

    mean(apply(cbind(low.excursion, up.excursion), 1, any))
  }

  # initialize the quantile piecewise linear function
  q0 = -Inf
  q  = rep(-Inf, length(x))
  I_weights = c(I_weights, Inf)

  # iterate over the different intervals
  for(k in 1:length(I_weights[-1])){
    if(I_weights[k] != 0){
      # find indices within the k-th interval
      subIk = subI[[k]]

      if(q0 == -Inf){
        # Get the local statistic
        maxIk = apply(rbind(apply(t(t(samples_minus[subIk,])), 2, max),
                            apply(t(t(samples_plus[subIk,])), 2, max)), 2, max)

        # Initialize the piecewise linear function
        mq = quantile( maxIk, 1 - alpha*I_weights[k], type = 8 )

        q[subIk] = mq

        # Get back into this loop if type is not linear
        if(type == "constant" || I_weights[k+1] == 0){
          q0 = -Inf
        }else{
          q0 = mq
        }
      }else{
        # define the function, which finds the optimal slope for the correct rejection rate
        solvef <- function(l){
          qq = q
          qq[subIk] = q0 - l * (x[subIk] - x[subIk[1]-1])

          return(subI_prob(k, qq) - alpha*I_weights[k])
        }
        qk <- uniroot(solvef, interval = c(-50, 50))

        # make the linear function only on the critical set
        ind = intersect(subIk, c(which(crit.set$minus), which(crit.set$plus)))

        q[ind] = q0 - qk$root * (x[ind] - x[subIk[1]-1])

        # Get back into this loop if type is not linear
        if(type == "constant" || I_weights[k+1] == 0){
          q0 = -Inf
        }else{
          q0 = q[subIk[length(subIk)]]
        }
      }
    }
  }

  I_weights = I_weights[-length(I_weights)]

  # Interval counter to later fill not important intervals
  if(any(is.infinite(q))){
    q[is.infinite(q)] = NA
  }

  # Get the linear function corresponding to the estimated values
  qq = approxfun(x, q, rule = 2, na.rm = TRUE)

  EmpRejections = IntervalProb(q = qq(x), crit.set = crit.set, samples = samples,
                               x = x, fair.intervals = knots, subI = subI)
  # plot(x, qq(x), type="l", main = "q function")
  # return the results
  return(list(u = qq, mu = mq, EmpRejections = EmpRejections))
}


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
fair_quantile_boot <- function(alpha, x, samples,
                               crit.set  = list(minus = rep(T, length(X)),
                                                plus  = rep(T, length(X))),
                               knots     = range(x),
                               I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                               type      = "linear",
                               alpha_up  = alpha*(length(knots) - 1),
                               maxIter = 20,
                               tol     = alpha / 10,
                               subI    = NULL,
                               inter   = NULL){
  # Get the interval
  if(is.null(subI) || is.null(inter)){
    s     = sub.intervals(x, knots, crit.set)
    subI  = s$subI
    inter = s$inter
  }

  # Spread the parts where no probability is expected fairly
  # to the other partitions
  if(!all(inter)){
    I_weights = I_weights + sum(I_weights[!inter]) / sum(inter)
    I_weights[!inter] = 0
  }

  # Initialize the u function
  alpha_k = alpha
  ufcns  <- fair_Bootstrap(alpha     = alpha_k, x = x, samples = samples,
                           crit.set  = crit.set,
                           knots     = knots,
                           I_weights = I_weights,
                           type      = type, subI = subI, inter = inter)

  diff <- ufcns$EmpRejections$global - alpha
  #  print(paste("Global rejection:", ufcns$EmpRejections$global))
  #  print(paste("alpha: ", alpha_k))

  niter   = 0
  if(maxIter != 0){
      a = c(alpha, alpha_up)

      while(niter < maxIter & abs(diff) > tol){
        if(niter < 3){
          alpha_k = a[1]*0.8 + a[2]*0.2
        }else{
          alpha_k = mean(a)
        }
        ufcns <- fair_Bootstrap(alpha     = alpha_k, x = x, samples = samples,
                                crit.set  = crit.set,
                                knots     = knots,
                                I_weights = I_weights,
                                type      = type, subI = subI, inter = inter)

        diff <- ufcns$EmpRejections$global - alpha

        if(diff < 0){
          a[1] = alpha_k
        }else{
          a[2] = alpha_k
        }
        niter = niter + 1
        #        print(paste("Global rejection:", ufcns$EmpRejections$global))
        #        print(paste("alpha: ", alpha_k))
      }
  }
#  print(paste("alpha QUIT: ", alpha_k))
  return(list(u = ufcns$u, du = ufcns$du, alpha_loc = alpha_k*I_weights,
              alpha_global = ufcns$EmpRejections$global, niter = niter))
}


#-------------------------------------------------------------------------------
# Kac Rice Algorithms
#-------------------------------------------------------------------------------
#' Find an optimal quantile function q in a basis by minimizing the L1/L2 norm
#' under constraints given by a Kac-Rice Formula.
#'
#' @param coef array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param basis an basis object from fda package
#' @return vector of length 2 containing the L1 norm and square of the L2 norm
#' of the function defined by coef and basis
#' @export
quantile_KRF <- function(x0, alpha, tau, type = "t", df = NULL, knots = c(0, 1),
                         alpha_weights = rep(1 / (length(knots) - 1), length(knots) - 1),
                         crossings = "up", lower.tail = FALSE, MAX = FALSE
                         ){
  # Initialize output
  summary <- list()

  # Get the necessary parameters from the input
  #-----------------------------------------------------------------------------
  # Get the amount of intervals
  NI = length(knots)-1
  # Get the correct KRF formula
  if(type == "t"){
    KRF <- function(u, du, x){
      KacRice_t( df  = df,
                 tau = tau,
                 u   = u,
                 du  = du,
                 x   = x,
                 crossings = "up", t0 = 1, sigma = 1,
                 lower.tail = FALSE, EC = TRUE)
    }
  }else if(type == "chi2"){
    KRF <- function(u, du, x){
      KacRice_chi2( df  = df,
                    tau = tau,
                    u   = u,
                    du  = du,
                    x   = x,
                    crossings = crossings, sigma = 1,
                    lower.tail = lower.tail, EC = TRUE)
    }
  }

  # Define the inital KRF constraint
  #-----------------------------------------------------------------------------
  KRF_constraint <- function(x){
    u0 <- x[1]
    m  <- x[-1]

    KRF(u  = function(x) eval_piecewise_u(u0, m, knots, x),
        du = function(x) eval_piecewise_du(u0, m, knots, x),
        x  = c(knots[1], knots[NI+1])) - alpha
  }

  ccc = ifelse(MAX, -1, 1) # change sign of target if max is needed

  # ----- Ziel + Gradient -----
  eval_f <- function(x) {
    list(
      objective = ccc*L2_cost(x, knots),           # deine J-Funktion (Skalar)
      gradient  = ccc*gr_L2_cost(x, knots)         # dein ∇J (oder numDeriv::grad(...))
    )
  }

  # ----- Gleichung g(x)=0 + numerischer Gradient -----
  eval_g_eq <- function(x) {
    g  <- as.numeric(KRF_constraint(x))          # Skalar!
    gg <- numDeriv::grad(function(z) as.numeric(KRF_constraint(z)), x)
    list(constraints = g, jacobian = matrix(gg, nrow = 1))
  }

  # ---- Gleichung g(x)=0 + schnelle Forward-Diff-Jacobian (1 x n) ----
  fast_grad_g <- function(x, eps_rel = 1e-6) {
    n  <- length(x)
    g0 <- as.numeric(KRF_constraint(x))
    gg <- numeric(n)
    for (i in 1:n) {
      h     <- eps_rel * (1 + abs(x[i]))
      xp    <- x; xp[i] <- xp[i] + h
      gi    <- as.numeric(KRF_constraint(xp))
      gg[i] <- (gi - g0) / h
    }
    gg
  }

  eval_g_eq <- function(x) {
    g  <- as.numeric(KRF_constraint(x))              # Skalar
    gg <- fast_grad_g(x)                              # Länge n
    list(constraints = g, jacobian = matrix(gg, nrow = 1))
  }

  # ----- Inequality constraint g(x)=0 + numerischer Gradient -----
  eval_g_ineq <- function(x) {
    eps = 0.5
    u0 <- x[1]
    m  <- x[-1]                     # Länge K
    dx <- diff(knots)               # Länge K

    # Werte von u an allen Knoten (K+1 Werte)
    u_at_knots <- c(u0, u0 + cumsum(m * dx))

    # g(x) <= 0 erzwingt u_at_knots >= eps
    g <- eps - u_at_knots

    # (A) Analytische Jacobimatrix (schnell & exakt):
    K <- length(dx)
    J <- matrix(0, nrow = K + 1, ncol = K + 1)  # (K+1 constraints) x (1+K Variablen)
    J[, 1] <- -1                                # d/du0 von (eps - u) = -1
    for (j in 1:K) {
      # m_j wirkt ab Knoten j+1 aufwärts mit Faktor -dx[j]
      J[(j + 1):(K + 1), j + 1] <- -dx[j]
    }

    list(constraints = g, jacobian = J)
  }

  # (optional) sinnvolle Box-Bounds zur Stabilisierung/Beschleunigung:
  K  <- length(knots) - 1

  lb <- x0[-1]
  lb[lb<0]  <- 1.3 * lb[lb<0]
  lb[lb>0]  <- 0.7 * lb[lb>0]
  lb[lb==0] <- 1.3 * mean(lb[lb<0])

  ub <- x0[-1]
  ub[ub<0]  <- 0.7 * ub[ub<0]
  ub[ub>0]  <- 1.3 * ub[ub>0]
  ub[ub==0] <- 1.3 * mean(ub[ub>0])

  lb = c(x0[1]*0.8, lb)
  ub = c(x0[1]*1.2, ub)

  if(all(diff(x0[-1])==0)){
    lb = c(x0[1]*0.8, -5/diff(knots))
    ub = c(x0[1]*1.2,  5/diff(knots))
  }

  if(type == "t"){
    res <- suppressWarnings(nloptr(
      x0 = x0,
      eval_f   = eval_f,
      eval_g_eq = eval_g_eq,
      lb = lb, ub = ub,
      opts = list(
        algorithm = "NLOPT_LD_SLSQP",
        maxeval   = 100,
        xtol_rel  = 1e-4,
        ftol_rel  = 1e-4,
        tol_constraints_eq = 1e-3,
        print_level = 0
      ))
    )
  }else{
    res <- suppressWarnings(nloptr(
      x0        = x0,
      eval_f    = eval_f,
      eval_g_eq = eval_g_eq,            # KRF-Gleichung
      eval_g_ineq = eval_g_ineq,        # Positivität an Knoten
      lb = lb, ub = ub,
      opts = list(
        algorithm = "NLOPT_LD_SLSQP",
        maxeval   = 100,               # ggf. etwas höher, da mehr Nebenbedingungen
        xtol_rel  = 1e-4,
        ftol_rel  = 1e-4,
        tol_constraints_eq   = 1e-3,
        tol_constraints_ineq = rep(1e-6, length(knots)),
        print_level = 0
      ))
    )
  }

  sqrt(L2_cost(res$solution, knots))
  KRF_constraint(res$solution)
  res$iterations

  compute_u_prefix(res$solution[1], res$solution[-1], knots)

  # Summarize results
  #-----------------------------------------------------------------------------
  u0 <- res$solution[1]
  m  <- res$solution[-1]

  u  = make_piecewise_u(u0, m, knots)
  du = make_piecewise_du(u0, m, knots)

  L2 <- suppressWarnings(sqrt(res$objective))
  L1 <- integrate_save(f1 = u, xlims = range(knots))
  width = c(L1, L2)
  names(width) <- c("L1", "L2")

  return(list(u = u, du = du, res = res,
              KRF_error = KRF_constraint(res$solution), width = width))
}

#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
fair_quantile_KRF <- function(alpha, tau, x = seq(0, 1, length.out = 2), df = NULL,
                              knots     = range(x),
                              I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                              u.type    = "linear",
                              type      = "t",
                              sigma     = 1,
                              alpha_up  = alpha*(length(knots)-1),
                              maxIter   = 20,
                              crossings = "up",
                              lower.tail = FALSE,
                              tol       = alpha / 10){
  # Get the correct KRF formula
  if(type == "t"){
    KRF <- function(u, du, x){
      KacRice_t( df = df,
                 tau = tau,
                 u   = u,
                 du  = du,
                 x   = x,
                 crossings = "up", t0 = 1, sigma = 1,
                 lower.tail = FALSE, EC = TRUE
      )
    }
  }else if(type == "chi2"){
    KRF <- function(u, du, x){
      KacRice_chi2( df  = df,
                    tau = tau,
                    u   = u,
                    du  = du,
                    x   = x,
                    crossings = crossings, sigma = 1,
                    lower.tail = lower.tail, EC = TRUE)
    }
  }
  # Get the correct algorithm 1 for the class of u and whether SCoPES or SCBs
  # are computed
  if(u.type == "linear" || u.type == "constant"){
    if(u.type == "linear"){
      alg1 <- function(alpha, tau){
        # alg1_KRF( alpha, tau = tau, df = df,
        #           type = type, knots = knots,
        #           I_weights = I_weights, sigma = sigma)
        alg1_fast( alpha, tau = tau, df = df,
                   type = type, knots = knots,
                   I_weights = I_weights, sigma = sigma,
                   crossings = crossings, lower.tail = lower.tail)
      }
    }else{
      alg1 <- function(alpha, tau){
        alg1_KRF_const( alpha, tau = tau, df = df,
                        type = type, knots = knots,
                        I_weights = I_weights, sigma = sigma)
      }
    }
    # Initialize the u function
    alpha_k = alpha
    ufcns  <- alg1(alpha_k, tau)

    # Find how far away u is from an alpha quantile
    diff <- KRF(u = ufcns$u, du = ufcns$du, x = range(knots)) - alpha

    niter   = 0
    if(abs(diff) > tol & maxIter != 0){
      alpha_k = alpha_up
      ufcns  <- alg1(alpha_k, tau)

      diff <- KRF(u = ufcns$u, du = ufcns$du, x = range(knots)) - alpha

      if(abs(diff) > tol){
        a = c(alpha, alpha_up)

        while(niter < maxIter & abs(diff) > tol){
          if(niter < 3){
            alpha_k = a[1]*0.8 + a[2]*0.2
          }else{
            alpha_k = mean(a)
          }

          ufcns <- alg1(alpha_k, tau)

          diff <- KRF(u = ufcns$u, du = ufcns$du, x = range(knots)) - alpha

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
  }else{
    KRF(u, du, x)
  }

  return(list(u = ufcns$u, du = ufcns$du, alpha_loc = alpha_k*I_weights,
              alpha_global = diff+alpha, niter = niter, width = ufcns$width ))
}


#' Find an optimal piecewise linear quantile function q to remove
#' conservativeness of standard Kac Rice formula approach for fair
#' thresholds.
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
alg1_KRF <- function(alpha, knots, tau, df = NULL, type = "t",
                     I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                     sigma = 1){
  # Get the correct KRF formula
  if(type == "t"){
    KRF <- function(u, du, x, crossings, t0){
      KacRice_t( df = df,
                 tau = tau,
                 u   = u,
                 du  = du,
                 x   = x,
                 crossings = crossings, t0 = t0, sigma = 1,
                 lower.tail = FALSE, EC = TRUE)
    }
  }
  # Get the amount of Intervals
  K <- length(I_weights)
  # Initialize the parameters for the piecewise linear function on each Interval
  u_fun <- matrix(NA, 2, K)

  ##############################################################################
  # Find constant for first interval
  find_u0 <- function(u0){
                  KRF(u  = function(y){rep(u0, length(y))},
                      du = function(y){rep(0, length(y))},
                      x  = c(knots[1], knots[2]),
                      crossings = "down", t0 = 2) - alpha * I_weights[1]
  }
  # slope =0 on the first interval
  u_fun[2, 1] <- 0
  # get the constant on the first interval
  u_fun[1, 1] <- suppressWarnings(uniroot(f = find_u0,
                                          interval = c(5e-2, 50),
                                          extendInt = "yes")$root)

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
      find_mk <- function(mk){
        KRF(u = function(y){u_fun[1, 2*k-2] + mk*(y - knots[2*k-2])},
            du = function(y){rep(mk, length(y))},
            x = c(knots[2*k-2], knots[2*k-1]),
            crossings = "up",
            t0 = 1) - alpha * I_weights[2*k-2]
      }
      u_fun[2, 2*k-2] <- suppressWarnings(uniroot(f = find_mk,
                                                  interval = c(-20, 20),
                                                  extendInt = "yes")$root)

      # String for p.linear function
      fct_body_u <- paste0(fct_body_u,
                           "+ c_v[",2*k-2,"]*pmin(pmax(t - knots[",2*k-2,"],0),  knots[",2*k-1,"] - knots[",2*k-2,"])")
      fct_body_du <- paste0(fct_body_du,
                            "+ c_v[",2*k-2,"]*pmax(sign(t - knots[",2*k-2,"]),0)*pmax(-sign(t - knots[",2*k-1,"]),0)")

      #---- odd intervals
      # Initialize first value on the interval
      u_fun[1, 2*k-1] <- u_fun[1, 2*k-2] +  u_fun[2, 2*k-2] * (knots[2*k-1] - knots[2*k-2])

      # Function to minimize to get slope
      find_mk2 <- function(mk){
        KRF(u  = function(y){u_fun[1, 2*k-1] + mk*(y - knots[2*k-1])},
            du = function(y){rep(mk, length(y))}, x = c(knots[2*k-1], knots[2*k]),
            crossings = "down", t0 = 2) - alpha * I_weights[2*k-1]
      }
      u_fun[2, 2*k-1] <- suppressWarnings(uniroot(f = find_mk2,
                                                  interval = c(-20, 20),
                                                  extendInt = "downX")$root)

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
    find_mK <- function(mK){
      KRF(u = function(y){u_fun[1, K] + mK*(y - knots[K])},
          du = function(y){rep(mK, length(y))},
          x = c(knots[K], knots[K+1]),
          crossings = "up", t0 = 1) - alpha * I_weights[K]
    }
    u_fun[2,K] <- suppressWarnings(uniroot(f = find_mK,
                                           interval = c(-20, 100),
                                           extendInt = "yes")$root)

    # String for p.linear function
    fct_body_u <- paste0(fct_body_u, "+ c_v[",K,"]*pmin(pmax(t - knots[",K,"],0),  knots[",K+1,"] - knots[",K,"])")
    fct_body_du <- paste0(fct_body_du, "+ c_v[",K,"]*pmax(sign(t - knots[",K,"]),0)*pmax(-sign(t - knots[",K+1,"]),0)")
  }

  coeffs = c(u_fun[1,1], u_fun[2,-1])

  ufun       <- function(t, c_v, knots){}
  body(ufun) <- parse(text = fct_body_u)
  dufun       <- function(t, c_v, knots){}
  body(dufun) <- parse(text = fct_body_du)

  u = function(t){ ufun(t, c_v = coeffs, knots = knots) }

  width <- c(integrate_save(f = Vectorize(function(x) abs(u(x)) ), c(0, 1) ),
             sqrt(L2_cost(c(u_fun[1,1], u_fun[2,]), knots)))
  names(width) <- c("L1", "L2")

  return(list(u = u, du = function(t){ dufun(t, c_v = u_fun[2,], knots = knots) },
              c_v = coeffs, width = width))
}

#' Find an optimal piecewise linear quantile function q to remove
#' conservativeness of standard Kac Rice formula approach for fair
#' thresholds.
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
alg1_fast <- function(alpha, knots, tau, df = NULL, type = "t",
                      I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                      crossings = "up", lower.tail = FALSE,
                      sigma = 1){


  # Get the correct KRF formula
  if(type == "t"){
    KRF <- function(u, du, x){
      KacRice_t( df  = df,
                 tau = tau,
                 u   = u,
                 du  = du,
                 x   = x,
                 crossings = "up", sigma = 1,
                 lower.tail = FALSE, EC = TRUE)
    }
  }else if(type == "chi2"){
    KRF <- function(u, du, x){
      KacRice_chi2( df  = df,
                 tau = tau,
                 u   = u,
                 du  = du,
                 x   = x,
                 crossings = crossings, sigma = 1,
                 lower.tail = lower.tail, EC = TRUE)
    }
  }


  # Get the amount of Intervals
  K <- length(I_weights)

  # Initialize the parameters for the piecewise linear function on each Interval
  u_fun <- matrix(NA, 2, K)

  ##############################################################################
  # Find constant for first interval
  find_u0 <- Vectorize(function(u0){
    KRF(u  = function(y){rep(u0, length(y))},
        du = function(y){rep(0, length(y))},
        x  = c(knots[1],  knots[2])) - alpha * I_weights[1]
  })
  # slope =0 on the first interval
  u_fun[2, 1] <- 0
  # get the constant on the first interval
  u_fun[1, 1] <- suppressWarnings(uniroot(f = find_u0,
                                          interval = c(1e-1, 50),
                                          extendInt = "yes")$root)
  #  u_fun[1, 1] <- bisect_from_negative(f = find_u0, interval = c(0, 50))$root

  # Initialize string for piecewise linear function
  fct_body_u  <- paste0("c_v[1]")
  fct_body_du <- paste0("c_v[1]")

  ##############################################################################
  # Find constant for first interval
  if(K == 1){
    width <- c(u_fun[1, 1], u_fun[1, 1])
    names(width) <- c("L1", "L2")

    return(list( u = Vectorize(function(t){u_fun[1, 1]}),
                 du = Vectorize(function(t){0}), width = width))
  }else if(K != 1){
    for(k in 2:K){
      #---- even intervals
      # Initialize start value on the interal
      u_fun[1, k] <- u_fun[1, k-1] + u_fun[2, k-1] * (knots[k] - knots[k-1])

      # Function to minimize to get slope
      find_mk <- Vectorize(function(mk){
        KRF(u = function(y){u_fun[1, k] + mk*(y - knots[k])},
            du = function(y){rep(mk, length(y))},
            x = c(knots[k], knots[k+1])) -  alpha * I_weights[k]
      })
      # Find the two options for the minimum
      u_fun[2, k] = suppressWarnings(uniroot(f = find_mk, interval = c(-20, 100),
                            tol = 1e-4, extendInt = "yes")$root)
      #     u_fun[2, k] = bisect_from_negative(f = find_mk,
      #                                         interval = c(-20, 20))$root

      # String for p.linear function
      fct_body_u <- paste0(fct_body_u,
                           "+ c_v[",k,"]*pmin(pmax(t - knots[",k,"],0),  knots[",k+1,"] - knots[",k,"])")
      fct_body_du <- paste0(fct_body_du,
                            "+ c_v[",k,"]*pmax(sign(t - knots[",k,"]),0)*pmax(-sign(t - knots[",k+1,"]),0)")
    }


    coeffs = c(u_fun[1,1], u_fun[2,-1])

    ufun       <- function(t, c_v, knots){}
    body(ufun) <- parse(text = fct_body_u)
    dufun       <- function(t, c_v, knots){}
    body(dufun) <- parse(text = fct_body_du)

    u = function(t){ ufun(t, c_v = coeffs, knots = knots) }

    width <- c(integrate_save(f = Vectorize(function(x) abs(u(x)) ), c(0, 1) ),
               sqrt(L2_cost(c(u_fun[1,1], u_fun[2,]), knots)))
    names(width) <- c("L1", "L2")

    return(list(u = u,
                du = function(t){
                  dufun(t, c_v = u_fun[2,], knots = knots) },
                c_v = coeffs, width = width))
  }
}

#' Find an optimal piecewise linear quantile function q to remove
#' conservativeness of standard Kac Rice formula approach for fair
#' thresholds.
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
quantile_nico <- function(alpha, knots, tau, df = NULL, type = "t",
                     I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                     weights = rep(1/(length(knots) - 1), length(knots) - 1),
                     sigma = 1){


  # Get the correct KRF formula
  if(type == "t"){
    KRF <- function(u, du, x, crossings, t0, EC){
      KacRice_t( df  = df,
                 tau = tau,
                 u   = u,
                 du  = du,
                 x   = x,
                 crossings = crossings, t0 = t0, sigma = 1,
                 lower.tail = FALSE, EC = EC)
    }

    Phi <- function(t){
                stats::pt(q = t, df = df, lower.tail = FALSE)
            }
  }
  # Get the amount of Intervals
  K <- length(I_weights)

  # Initialize the parameters for the piecewise linear function on each Interval
  u_fun <- matrix(NA, 2, K)

  ##############################################################################
  # Find constant for first interval
  find_u0 <- function(u0){
    KRF(u  = function(y){rep(u0, length(y))},
        du = function(y){rep(0, length(y))},
        x  = c(knots[1],  knots[2]),
        crossings = "up", t0 = 1, EC = F) + Phi(u0) * weights[1] - alpha * I_weights[1]
  }
  # slope =0 on the first interval
  u_fun[2, 1] <- 0
  # get the constant on the first interval
  u_fun[1, 1] <- suppressWarnings(uniroot(f = find_u0,
                                          interval = c(5e-2, 50),
                                          extendInt = "upX")$root)
#  u_fun[1, 1] <- bisect_from_negative(f = find_u0, interval = c(0, 50))$root

  # Initialize string for piecewise linear function
  fct_body_u  <- paste0("c_v[1]")
  fct_body_du <- paste0("c_v[1]")

  ##############################################################################
  # Find constant for first interval
  if(K == 1){
    width <- c(u_fun[1, 1], u_fun[1, 1])
    names(width) <- c("L1", "L2")

    return(list( u = Vectorize(function(t){u_fun[1, 1]}),
                du = Vectorize(function(t){0}), width = width))
  }else if(K != 1){
    for(k in 2:K){
      #---- even intervals
      # Initialize start value on the interal
      u_fun[1, k] <- u_fun[1, k-1] + u_fun[2, k-1] * (knots[k] - knots[k-1])

      # Function to minimize to get slope
      find_mk <- Vectorize(function(mk){
        KRF(u = function(y){u_fun[1, k] + mk*(y - knots[k])},
            du = function(y){rep(mk, length(y))},
            x = c(knots[k], knots[k+1]),
            crossings = "up",
            t0 = NULL, EC = FALSE)  + Phi(u_fun[1, 1])*weights[k] -
                alpha * I_weights[k]
      })
      # Find the two options for the minimum
      u_fun[2, k] = suppressWarnings(uniroot(f = find_mk,
                                             interval = c(-20, 100),
                                             extendInt = "yes", tol = 1e-4)$root)
 #     u_fun[2, k] = bisect_from_negative(f = find_mk,
#                                         interval = c(-20, 20))$root

      # String for p.linear function
      fct_body_u <- paste0(fct_body_u,
                           "+ c_v[",k,"]*pmin(pmax(t - knots[",k,"],0),  knots[",k+1,"] - knots[",k,"])")
      fct_body_du <- paste0(fct_body_du,
                            "+ c_v[",k,"]*pmax(sign(t - knots[",k,"]),0)*pmax(-sign(t - knots[",k+1,"]),0)")
    }


  coeffs = c(u_fun[1,1], u_fun[2,-1])

  ufun       <- function(t, c_v, knots){}
  body(ufun) <- parse(text = fct_body_u)
  dufun       <- function(t, c_v, knots){}
  body(dufun) <- parse(text = fct_body_du)

  u = function(t){ ufun(t, c_v = coeffs, knots = knots) }

  width <- c(integrate_save(f = Vectorize(function(x) abs(u(x)) ), c(0, 1) ),
             sqrt(L2_cost(c(u_fun[1,1], u_fun[2,]), knots)))
  names(width) <- c("L1", "L2")

  # local bounds
  loc_bds <- Phi(u(knots[-length(knots)])) - weights * Phi(u(knots[1]))
  print(alpha*knots[2] + colMeans(matrix(loc_bds, nrow = 8 )))

  return(list(u = u,
              du = function(t){
              dufun(t, c_v = u_fun[2,], knots = knots) },
              c_v = coeffs, width = width,
              loc_bds = loc_bds))
  }
}


#' Find an optimal piecewise constant quantile function q to remove
#' conservativeness of standard Kac Rice formula approaches.
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
alg1_KRF_const <- function(alpha, tau, df = NULL, type = "t", knots,
                           I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                           sigma = 1){
  # Get the correct KRF formula
  if(type == "t"){
    KRF <- function(u, du, x, crossings, t0){
      KacRice_t( df = df,
                 tau = tau,
                 u   = u,
                 du  = du,
                 x   = x,
                 crossings = crossings, t0 = t0, sigma = 1,
                 lower.tail = FALSE, EC = TRUE
      )
    }
  }
  # Get the amount of Intervals
  K <- length(I_weights)
  # Initialize the parameters for the piecewise linear function on each Interval
  u_fun <- rep(NA, K)
  fct_body_u  <- paste0("0")
  fct_body_du <- paste0("0")

  ##############################################################################
  # Find constant for first interval
  for(k in 1:K){
    find_u0 <- function(u0){
      KRF(u = function(y){rep(u0, length(y))},
          du = function(y){rep(0, length(y))},
          x = c(knots[k], knots[k+1]),
          crossings = "up", t0 = 1) - alpha * I_weights[k]
    }
    # get the constant on the first interval
    u_fun[k] <- suppressWarnings(uniroot(f = find_u0,
                                         interval = c(1e-1, 50),
                                         extendInt = "downX")$root)

    # Initialize string for piecewise linear function
    fct_body_u <- paste0(fct_body_u, "+ c_v[",k,"]*ifelse(sign(t - knots[",k,"])==0, 1, pmax(sign(t - knots[",k,"]),0))*pmax(-sign(t - knots[",k+1,"]),0)")
    fct_body_du <- paste0(fct_body_du, "+ c_v[",k,"]")
  }
  fct_body_u <- paste0(fct_body_u, "+ c_v[",K,"]*ifelse(sign(t - knots[",K+1,"])==0, 1, 0)")
  ufun        <- function(t, c_v, knots){}
  body(ufun)  <- parse(text = fct_body_u)
  dufun       <- function(t, c_v, knots){}
  body(dufun) <- parse(text = fct_body_du)

  return(list( u = function(t){  ufun(t, c_v = u_fun, knots = knots) },
               du = function(t){ dufun(t, c_v = rep(0,K), knots = knots) }))
}


#' This function is only here for comparison with the implementation of
#' D Liebl. It is essentially copied from the ffscb package.
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
alg1_z_DL <- function (tau, diag.cov=rep(1, length(tau)), conf.level = 0.95, n_int = 4){
  tol <- .Machine$double.eps^0.35
  alpha <- 1 - conf.level
  tt <- seq(0, 1, len = length(tau))
  tau_v <- tau
  tau_f <- function(t) {
    stats::approx(x = seq(0, 1, len = length(tau)), y = tau,
                  xout = t)$y
  }
  knots <- seq(0, 1, len = (n_int + 1))
  if (!is.numeric(n_int)) {
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  if (n_int <= 0) {
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  if (n_int%%1 != 0) {
    stop("n_int must be a stricitly positive integer value (n_int=1,2,3,...)")
  }
  if (n_int == 1) {
    tau01 <- sum(tau_v) * diff(tt)[1]
    myfun1 <- function(c1) {
      stats::pnorm(q = c1, lower.tail = F) + exp(-c1^2/2) *
        tau01/(2 * pi) - (alpha/2)
    }
    const_band <- stats::uniroot(f = myfun1, interval = c(0,
                                                          10), extendInt = "downX")$root
    band <- const_band * sqrt(diag.cov)
    return(band)
  }
  c_v <- numeric(n_int)
  const_int <- 1
  fct_body <- paste0("c_v[", const_int, "]")
  for (j in (const_int + 1):n_int) {
    fct_body <- paste0(fct_body, "+ c_v[", j, "]*pmax(t - knots[",
                       j, "],0)")
  }
  ufun <- function(t, c_v, knots) {
  }
  body(ufun) <- parse(text = fct_body)
  tau_init <- sum(tau_v[knots[const_int] <= tt & tt <= knots[const_int +
                                                               1]]) * diff(tt)[1]
  myfun1 <- function(c1) {
    stats::pnorm(-c1) + c(exp(-c1^2/2) * tau_init/(2 * pi) -
                            (alpha/2)/n_int)
  }
  c_v[const_int] <- stats::uniroot(f = myfun1, interval = c(0,
                                                            10), extendInt = "downX", tol = tol)$root
  for (j in (const_int + 1):n_int) {
    myfunj <- function(cj) {
      if (j == (const_int + 1)) {
        c_v_sum <- 0
      }
      else {
        c_v_sum <- c_v[(const_int + 1):(j - 1)]
      }
      ufun_j <- function(t, cj) {
        ufun(t = t, c_v = c(c_v[const_int:(j - 1)], cj,
                            rep(0, times = (n_int - j))), knots = knots)
      }
      fn1 <- function(t, cj) {
        (tau_f(t)/(2 * pi)) * exp(-ufun_j(t, cj)^2/2) *
          exp(-sum(c(c_v_sum, cj))^2/(2 * tau_f(t)^2))
      }
      fn2 <- function(t, cj) {
        sum(c(c_v_sum, cj))/sqrt(2 * pi) * exp(-ufun_j(t,
                                                       cj)^2/2) * stats::pnorm(sum(c(c_v_sum, cj))/tau_f(t))
      }
      fn3 <- function(t, cj) {
        sum(c(c_v_sum, cj))/sqrt(2 * pi) * exp(-ufun_j(t,
                                                       cj)^2/2) * stats::pnorm(-sum(c(c_v_sum, cj))/tau_f(t))
      }
      intgr1 <- sum(fn1(t = tt[knots[j] < tt & tt <= knots[j +
                                                             1]], cj = cj)) * diff(tt)[1]
      if (j%%2 != 0) {
        intgr2 <- sum(fn2(t = tt[knots[j] < tt & tt <=
                                   knots[j + 1]], cj = cj)) * diff(tt)[1]
      }
      else {
        intgr2 <- 0
      }
      if (j%%2 == 0) {
        intgr3 <- sum(fn3(t = tt[knots[j] < tt & tt <=
                                   knots[j + 1]], cj = cj)) * diff(tt)[1]
      }
      else {
        intgr3 <- 0
      }
      if (j%%2 == 0) {
        res <- c(stats::pnorm(-ufun(knots[j], c_v = c_v,
                                    knots = knots)) + intgr1 + intgr2 - intgr3 -
                   (alpha/2)/n_int)
      }
      if (j%%2 != 0) {
        res <- c(stats::pnorm(-ufun_j(t = knots[j + 1],
                                      cj)) + intgr1 + intgr2 - intgr3 - (alpha/2)/n_int)
      }
      return(res)
    }
    c_v[j] <- stats::uniroot(f = myfunj, interval = c(-10,
                                                      10), extendInt = "downX", tol = tol)$root
  }
  band.eval <- ufun(t = tt, c_v = c_v, knots = knots)
  band <- band.eval * sqrt(diag.cov)
  return(band)
}

