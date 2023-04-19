#------------------------------------------------------------------------------#
#                                                                              #
#     Bootstrap functions to obtain the quantiles of the maximum of random     #
#     fields                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#      - Fair_quantile_boot()
#      - Optimize_Fair_quantile_boot()
#      - fair_Bootstrap()
#      - fair_quantile_boot()
#      - alg1_gKR_t()
#      - alg1_gKR_const()
#      - fair_quantile_EEC_t()
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
Fair_quantile_boot <- function(samples,
                               x,
                               fair.intervals,
                               fair.type = "linear",
                               crit.set  = list(minus = rep(T, length(x)),
                                                plus  = rep(T, length(x))),
                               alpha     = 0.05,
                               diff.fair = NULL,
                               subI      = NULL,
                               inter     = NULL ){
  #
  if(is.null(subI) || is.null(inter)){
    s = sub.intervals(x, fair.intervals, crit.set)
    subI  = s$subI
    inter = s$inter
  }

  #
  Sminus = crit.set$minus
  Splus  = crit.set$plus

  # Get the length of the different intervals for the fairness
  if(is.null(diff.fair)){
    dfair = diff(fair.intervals)
  }else{
    dfair = diff.fair
  }

  # Spread the parts where no probability is expected fairly
  # to the other partitions
  if(!all(inter)){
    dfair = dfair + sum(dfair[!inter]) / sum(inter)
    dfair[!inter] = 0
  }

  dimS = dim(samples)
  #  if(fair.type = "linear"){
  samples_minus = -samples
  samples_minus[!Sminus,] <- -Inf
  samples_plus = samples
  samples_plus[!Splus,]   <- -Inf
  #}

  subI_prob <- function(k, q){
    if(is.null(crit.set$minus)){
      low.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_minus     = intersect(which(crit.set$minus), subI[[k]])
      if(length(ind_minus) == 0){
        low.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_minus) == 1){
        low.excursion = -samples[ind_minus,] - q[ind_minus] >= 0
      }else{
        low.excursion = apply(-samples[ind_minus,] - q[ind_minus], 2, max) >= 0
      }
    }

    if(is.null(crit.set$plus)){
      up.excursion = rep(FALSE, dim(samples)[2])
    }else{
      ind_plus     = intersect(which(crit.set$plus), subI[[k]])
      if(length(ind_plus) == 0){
        up.excursion = rep(FALSE, dim(samples)[2])
      }else if(length(ind_plus) == 1){
        up.excursion = samples[ind_plus,] - q[ind_plus] >= 0
      }else{
        up.excursion = apply(samples[ind_plus,] - q[ind_plus], 2, max) >= 0      }
    }

    mean(apply(cbind(low.excursion, up.excursion), 1, any))
  }

  # initialize the quantile piecewise linear function
  q0 = -Inf
  q  = rep(-Inf, length(x))

  # Get the indices of the partitions which are used
  # to control the limit distribution
  Ik = c(which(inter), -Inf)

  # iterate over the different intervals
  for(k in 1:length(Ik[-length(Ik)])){
    # find indices within the k-th interval
    subIk = subI[[Ik[k]]]

    if(q0 == -Inf){
      # Get the local test statistic
      maxIk = apply(rbind(apply(samples_minus[subIk,], 2, max),
                          apply(samples_plus[subIk,], 2, max)), 2, max)

      # Initialize the piecewise linear function
      mq = quantile( maxIk, 1-alpha*dfair[Ik[k]], type = 8 )

      q[subIk] = mq

      # Get back into this loop if type is not linear
      if(fair.type == "constant" || Ik[k]+1 != Ik[k+1]){
        q0 = -Inf
      }else{
        q0 = mq
      }
    }else{
      # define the function, which finds the optimal slope for the correct rejection rate
      solvef <- function(l){
        qq = q
        qq[subIk] = q0 - l * (x[subIk] - x[subIk[1]-1])

        return(subI_prob(Ik[k], qq) - alpha * dfair[Ik[k]])
      }
      qk <- uniroot(solvef, interval = c(-500, 500))

      q[subIk] = q0 - qk$root * (x[subIk] - x[subIk[1]-1])

      # Get back into this loop if type is not linear
      if(fair.type == "constant" || Ik[k]+1 != Ik[k+1]){
        q0 = -Inf
      }else{
        q0 = q[subIk[length(subIk)]]
      }
    }
  }

  # Interval counter to later fill not important intervals
  if(any(is.infinite(q))){
    if(fair.type == "linear"){
      q[is.infinite(q)] = NA
      if(is.na(q[1])){
        q[1] = mean(q, na.rm = TRUE)
      }
      Eqx = length(q)
      if(is.na(q[Eqx])){
        q[Eqx] = mean(q, na.rm = TRUE)
      }
      q = na.approx(q)
    }else{
      q[is.infinite(q)] = mean(q[!is.infinite(q)])
    }
  }

  EmpRejections = IntervalProb(q, crit.set, samples, x, fair.intervals, subI = subI)

  # return the results
  return(list(q = q, mq = mq, EmpRejections = EmpRejections))
}


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
Optimize_Fair_quantile_boot <- function(samples,
                                    x,
                                    fair.intervals,
                                    fair.type = "linear",
                                    crit.set,
                                    alpha  = 0.05,
                                    niter  = 10,
                                    diff.fair = NULL,
                                    subI  = NULL,
                                    inter = NULL,
                                    print.coverage = TRUE ){

  # Fill subI and inter, if not provided
  if(is.null(subI) || is.null(inter)){
    s = sub.intervals(x, fair.intervals, crit.set)
    subI  = s$subI
    inter = s$inter
  }

  # Get the stepsize for the interval iteration
  if( sum(crit.set$minus)/2 + sum(crit.set$minus)/2 > 0.75*length(x)  ){
    stepsize = 3
  }else{
    stepsize = 1.5
  }

  # Compute the fair threshold function
  test = Fair_quantile_boot(samples = samples,
                         x = x,
                         fair.intervals = fair.intervals,
                         fair.type = fair.type,
                         crit.set  = crit.set,
                         alpha     = alpha,
                         diff.fair = diff.fair,
                         subI = subI, inter = inter)

  # Initialize values for computing the fair threshold function
  count = 0
  breakCond = FALSE
  oldEmp    = test$EmpRejections$global
  eps = max(alpha * 0.05, 10/dim(samples)[2])
  alpha_new = alpha

  # Loop to remove conservativeness of the fair threshold function
  # using nested intervals
  while(breakCond == FALSE && count < niter){
    diffCoverage = oldEmp - alpha
    if( abs(diffCoverage) > eps ){
      if(count == 0){
        if(diffCoverage < 0){
          a = c(alpha, stepsize * alpha)
        }else{
          a = c(0.1 * alpha, alpha)
        }
      }else{
        if(diffCoverage < 0){
          a[1] = alpha_new
          if(a[1] == a[2]){
            a[2] = 2 * a[2]
          }
        }else{
          a[2] = alpha_new
          if(a[1] == a[2]){
            a[1] = 0.1 * a[1]
          }
        }
      }

      # update the global alpha
      alpha_new = mean(a)

      # Get new quantile function
      test   = Fair_quantile_boot(samples = samples,
                                  x = x,
                                  fair.intervals = fair.intervals,
                                  fair.type = fair.type,
                                  crit.set  = crit.set,
                                  alpha     = alpha_new,
                                  diff.fair = diff.fair,
                                  subI  = subI,
                                  inter = inter )
      count  = count + 1
      oldEmp = test$EmpRejections$global
      if(print.coverage){
        print(oldEmp)
      }
    }else{
      breakCond = TRUE
    }
  }

  # return the results
  return(test)
}


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
fair_Bootstrap <- function(alpha,
                           samples,
                           x,
                           knots,
                           type = "linear",
                           crit.set  = list(minus = rep(T, length(x)),
                                            plus  = rep(T, length(x))),
                           I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
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
        up.excursion = apply(samples[ind_plus,] - q[ind_plus], 2, max) >= 0      }
    }

    mean(apply(cbind(low.excursion, up.excursion), 1, any))
  }

  # initialize the quantile piecewise linear function
  q0 = -Inf
  q  = rep(-Inf, length(x))

  # Get the indices of the partitions which are used
  # to control the limit distribution
  Ik = c(which(inter), -Inf)

  # iterate over the different intervals
  for(k in 1:length(Ik[-length(Ik)])){
    # find indices within the k-th interval
    subIk = subI[[Ik[k]]]

    if(q0 == -Inf){
      # Get the local test statistic
      maxIk = apply(rbind(apply(samples_minus[subIk,], 2, max),
                          apply(samples_plus[subIk,], 2, max)), 2, max)

      # Initialize the piecewise linear function
      mq = quantile( maxIk, 1 - alpha*I_weights[Ik[k]], type = 8 )

      q[subIk] = mq

      # Get back into this loop if type is not linear
      if(type == "constant" || Ik[k]+1 != Ik[k+1]){
        q0 = -Inf
      }else{
        q0 = mq
      }
    }else{
      # define the function, which finds the optimal slope for the correct rejection rate
      solvef <- function(l){
        qq = q
        qq[subIk] = q0 - l * (x[subIk] - x[subIk[1]-1])

        return(subI_prob(Ik[k], qq) - alpha*I_weights[Ik[k]])
      }
      qk <- uniroot(solvef, interval = c(-500, 500))

      q[subIk] = q0 - qk$root * (x[subIk] - x[subIk[1]-1])

      # Get back into this loop if type is not linear
      if(type == "constant" || Ik[k]+1 != Ik[k+1]){
        q0 = -Inf
      }else{
        q0 = q[subIk[length(subIk)]]
      }
    }
  }

  # Interval counter to later fill not important intervals
  if(any(is.infinite(q))){
    q[is.infinite(q)] = NA
  }

  # Get the linear function corresponding to the estimated values
  qq = approxfun(x, q, rule = 2, na.rm = TRUE)

  EmpRejections = IntervalProb(qq(x), crit.set, samples, x, knots, subI = subI)

  # return the results
  return(list(u = qq, mu = mq, EmpRejections = EmpRejections))
}


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
fair_quantile_boot <- function(alpha, df = NULL, knots, samples, sigma = 1,
                                I_weights = rep(1/(length(knots) - 1), length(knots) - 1),
                                alpha_up  = alpha*(length(knots) - 1),
                                maxIter = 20,
                                tol     = alpha / 100,
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
  ufcns  <- fair_Bootstrap(alpha = alpha_k, df = df, knots = knots,
                           samples = samples, sigma = sigma,
                           I_weights = I_weights,
                           subI = subI, inter = inter)

  diff <- ufcns$EmpRejections$global - alpha

  niter   = 0
  if(abs(diff) > tol & maxIter != 0){

    alpha_k = alpha_up
    ufcns  <- fair_Bootstrap(alpha = alpha_k, df = df, knots = knots,
                             samples = samples, sigma = sigma,
                             I_weights = I_weights,
                             subI = subI, inter = inter)

    diff <- ufcns$EmpRejections$global - alpha

    if(abs(diff) > tol){
      a = c(alpha, alpha_up)

      while(niter < maxIter & abs(diff) > tol){
        alpha_k = a[1]*0.6 + a[2]*0.4
        ufcns <- fair_Bootstrap(alpha = alpha_k, df = df, knots = knots,
                                samples = samples, sigma = sigma,
                                I_weights = I_weights,
                                subI = subI, inter = inter)

        diff <- ufcns$EmpRejections$global - alpha

        if(diff < 0){
          a[1] = alpha_k
        }else{
          a[2] = alpha_k
        }
        niter = niter + 1
      }
    }
  }

  return(list(u = ufcns$u, du = ufcns$du, alpha_loc = alpha_k*I_weights, niter = niter))
}


#-------------------------------------------------------------------------------
# Kac Rice Algorithms
#-------------------------------------------------------------------------------
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
alg1_gKR_t <- function(alpha, knots, tau, df = NULL, sigma = 1,
                       I_weights = rep(1/(length(knots) - 1), length(knots) - 1)){
  # Get the amount of Intervals
  K <- length(I_weights)
  # Initialize the parameters for the piecewise linear function on each Interval
  u_fun <- matrix(NA, 2, K)

  ##############################################################################
  # Find constant for first interval
  find_u0 <- function(u0){
    GeneralKacRice_t(tau = tau,
          u = function(y){rep(u0, length(y))},
          du = function(y){rep(0, length(y))}, x = c(knots[1], knots[2]),
          df = df,
          sigma = sigma, crossings = "down", t0 = 2) - alpha * I_weights[1]
  }
  # slope =0 on the first interval
  u_fun[2, 1] <- 0
  # get the constant on the first interval
  u_fun[1, 1] <- uniroot(f = find_u0, interval = c(0, 50), extendInt = "downX")$root

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
        GeneralKacRice_t(tau, u = function(y){u_fun[1, 2*k-2] + mk*(y - knots[2*k-2])},
              du = function(y){rep(mk, length(y))}, x = c(knots[2*k-2], knots[2*k-1]),
              sigma = sigma, df = df, crossings = "up", t0 = 1) - alpha * I_weights[2*k-2]
      }
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
      find_mk2 <- function(mk){
        GeneralKacRice_t(tau,
              u  = function(y){u_fun[1, 2*k-1] + mk*(y - knots[2*k-1])},
              du = function(y){rep(mk, length(y))}, x = c(knots[2*k-1], knots[2*k]),
              sigma = sigma, df = df, crossings = "down", t0=2) - alpha * I_weights[2*k-1]
      }
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
    find_mK <- function(mK){
      GeneralKacRice_t(tau, u = function(y){u_fun[1, K] + mK*(y - knots[K])},
            du = function(y){rep(mK, length(y))}, x = c(knots[K], knots[K+1]),
            df = df, sigma = sigma, crossings = "up") - alpha * I_weights[K]
    }
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
alg1_gKR_const <- function(alpha, knots, tau, df = NULL, sigma = 1,
                         I_weights = rep(1/(length(knots) - 1), length(knots) - 1)){
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
      GeneralKacRice_t(tau = tau,
            u = function(y){rep(u0, length(y))},
            du = function(y){rep(0, length(y))}, x = c(knots[k], knots[k+1]),
            df = df, sigma = sigma, crossings = "up", t0 = 1) - alpha * I_weights[k]
    }
    # get the constant on the first interval
    u_fun[k] <- uniroot(f = find_u0, interval = c(0, 50), extendInt = "downX")$root

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
  ufcns  <- alg1_gKR_t(alpha = alpha_k, df = df, knots = knots,
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
    ufcns  <- alg1_gKR_t(alpha = alpha_k, df = df, knots = knots,
                         tau = tau, sigma = sigma, I_weights = I_weights)

    diff <- GeneralKacRice_t(tau = tau, u = ufcns$u, du = ufcns$du,
                             df = df, x = range(knots),
                             crossings = "up", t0 = 1) - alpha

    if(abs(diff) > tol){
      a = c(alpha, alpha_up)

      while(niter < maxIter & abs(diff) > tol){
        alpha_k = a[1]*0.6 + a[2]*0.4
        ufcns <- alg1_gKR_t(alpha = alpha_k, df = df, knots = knots,
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
alg1_z_DL <- function (tau, diag.cov, conf.level = 0.95, n_int = 4){
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
