#------------------------------------------------------------------------------#
#                                                                              #
#     Auxillary functions for the package
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
# - IntervalProb()
# - sub.intervals()
# - q_method_gen
# - Preimage_method_gen()
# - plot_col()
#------------------------------------------------------------------------------#
# Developer notes:
# - Fix the documentation and man pages
#
#------------------------------------------------------------------------------#
#' Estimates using a sample of random functions (for example a bootstrap sample)
#' the FWER on different intervals.
#'
#' @param sample array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param x vector of length K of locations at which the sample is observed.
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
IntervalProb <- function(q, crit.set, samples, x, fair.intervals, subI = NULL){
  # Get the subinterval indices
  if(is.null(subI)){
    subI = list()
    for(k in 2:length(fair.intervals)){
      subI[[k-1]] = which(x  >= fair.intervals[k-1] & x  <= fair.intervals[k])
    }
  }

  subI_prob <- function(k){
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
        up.excursion = apply(samples[ind_plus,] - q[ind_plus], 2, max) >= 0
      }
    }

    mean(apply(cbind(low.excursion, up.excursion), 1, any))
  }

  # Intervalwise rejections
  subI.probs = vapply(1:length(subI), subI_prob, FUN.VAL = 0.1)

  # Global rejectionrate
  if(is.null(crit.set$minus)){
    low.excursion = rep(FALSE, dim(samples)[2])
  }else{
    if(sum(crit.set$minus) == 0){
      low.excursion = rep(FALSE, dim(samples)[2])
    }else if(sum(crit.set$minus) == 1){
      low.excursion = -samples[crit.set$minus,] - q >= 0
    }else{
      low.excursion = apply(-samples[crit.set$minus,] - q, 2, max) >= 0
    }

  }
  if(is.null(crit.set$plus)){
    up.excursion = rep(FALSE, dim(samples)[2])
  }else{
    if(sum(crit.set$plus) == 0){
      up.excursion =  rep(FALSE, dim(samples)[2])
    }else if(sum(crit.set$plus) == 1){
      up.excursion = samples[crit.set$plus,] - q >= 0
    }else{
      up.excursion = apply(samples[crit.set$plus,] - q, 2, max) >= 0
    }
  }
  global.prob = mean(apply(cbind(low.excursion, up.excursion), 1, any))

  return(list(global = global.prob, local = subI.probs))
}

#' Estimates using a sample of random functions (for example a bootstrap sample)
#' the FWER on different intervals.
#'
#' @param sample array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param x vector of length K of locations at which the sample is observed.
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
sub.intervals <- function(x, fair.intervals, crit.set){
  # Get the subintervals
  subI = list()
  for(k in 2:length(fair.intervals)){
    subI[[k-1]] = which(x >= fair.intervals[k-1] & x <= fair.intervals[k])
  }

  # Get the intervals which have an intersection with the critical sets
  subI.cit <- unlist(lapply(subI, function(l){
    ifelse(length(intersect(l, which(crit.set$minus))) != 0 |
             length(intersect(l, which(crit.set$plus))) != 0, TRUE, FALSE )

  } ))

  list(subI = subI, inter = subI.cit)
}

#' This function generates a list with the required input for different
#' methods to estimate the quantile function q of a SCoPE set.
#'
#' @param name String name of the implemented method to
#' estimate the quantile function
#' @param ... Optional parameters for the quantile estimation method
#' @return q.method
#' @export
q_method_gen <- function(name, ...){
  # Get the input arguments
  arguments = as.list(match.call())
  arguments = names(arguments)

  # Initialize the q_method list, which will be the output
  q.method <- list(
    name = name
  )

  # Fill the q_method list for different quantile function estimation
  # methods
  if(q.method$name %in% c("mboot", "fair.mboot")){
    # Use default values for the multiplier bootstrap, if
    # the values are not provided in the input arguments.
    if(!("Boottype" %in% arguments)){
      Boottype = "t"
    }
    if(!("weights" %in% arguments)){
      weights = "rademacher"
    }
    if(!("Mboots" %in% arguments)){
      Mboots = 5e3
    }

    q.method$Boottype   = Boottype
    q.method$weights    = weights
    q.method$Mboots     = Mboots
    q.method$R          = R

    # Add the fairness parameters if the bootstrap
    # should be fair
    if(q.method$name == "fair.mboot"){
      q.method$fair.intervals = fair.intervals
      if(!("fair.type" %in% arguments)){
        q.method$fair.type = "linear"
      }else{
        q.method$fair.type = fair.type
      }

      if(!("fair.niter" %in% arguments)){
        q.method$fair.niter = 10
      }else{
        q.method$fair.niter = fair.niter
      }
    }
  }else if(q.method$name == "gKR_t"){
    if(!("df" %in% arguments)){
      q.method$df = NULL
    }else{
      q.method$df = df
    }
    if(!("knots" %in% arguments)){
      q.method$knots = c(0,1)
    }else{
      q.method$knots = knots
    }
    q.method$Nknots = length(q.method$knots) - 1
    if(!("tau.est" %in% arguments)){
      q.method$tau.est = tau_est
    }else{
      q.method$tau.est = tau.est
    }
    if(!("I_weights" %in% arguments)){
      q.method$I_weights = rep(1/(length(q.method$knots) - 1), length(q.method$knots) - 1)
    }else{
      q.method$I_weights = I_weights
    }
    if(!("alpha_up" %in% arguments)){
      q.method$alpha_up = alpha*(length(q.method$knots)-1)
    }else{
      q.method$alpha_up = alpha_up
    }
    if(!("maxIter" %in% arguments)){
      q.method$maxIter = 0
    }else{
      q.method$maxIter = maxIter
    }

  }else if(q.method$name == "t.iid"){
    q.method$df      = N - 1
  }

  if("print.coverage" %in% arguments){
    q.method$print.coverage = print.coverage
  }else{
    q.method$print.coverage = TRUE
  }

  return(q.method)
}

#' This function generates a list with the required input for different
#' methods to estimate the quantile function q of a SCoPE set.
#'
#' @param name String name of the implemented method to
#' estimate the quantile function
#' @param ... Optional parameters for the quantile estimation method
#' @return Scale field
#' @export
Preimage_method_gen <- function(name, ...){
  # Get the input arguments
  #  arguments = as.list(match.call())

  # Initialize the Preimage.method list, which will be the output
  Preimage.method <- list(
    name = name
  )

  # Fill the Preimage.method list for different Preimage estimation
  # methods
  if(Preimage.method$name == "thickening"){
    Preimage.method$kN = kN

  }else if(Preimage.method$name == "true"){
    Preimage.method$kN = 0
    Preimage.method$mu = mu
  }else if(Preimage.method$name == "SCB"){
    Preimage.method$kN = 0
  }else if(Preimage.method$name == "Storey.iid"){
    Preimage.method$m0 = m0
  }

  return(Preimage.method)
}


#' @export
plot_col <- function(statistic, C, detect, x = NULL,
                     xlab = '', ylab = '', title = '',
                     mu = NULL){
  if(is.null(x)){
    x = 1:length(statistic)
  }
  y = statistic
  # Get a color vector indicating the "red" (upper excursions) and
  # the blue (lower excursions) set
  colVec <- rep("black", length(y))
  colVec[detect] <- "red"
    plot(x, y, col = colVec,
         pch = 18, xlab = xlab, ylab = ylab, main = title)
    lines(x, C, lty = 2, col = "orchid3")
}

