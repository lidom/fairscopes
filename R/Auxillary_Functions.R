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
      low.excursion = -samples[crit.set$minus,] - q[crit.set$minus] >= 0
    }else{
      low.excursion = apply(-samples[crit.set$minus,] - q[crit.set$minus], 2, max) >= 0
    }

  }
  if(is.null(crit.set$plus)){
    up.excursion = rep(FALSE, dim(samples)[2])
  }else{
    if(sum(crit.set$plus) == 0){
      up.excursion =  rep(FALSE, dim(samples)[2])
    }else if(sum(crit.set$plus) == 1){
      up.excursion = samples[crit.set$plus,] - q[crit.set$plus] >= 0
    }else{
      up.excursion = apply(samples[crit.set$plus,] - q[crit.set$plus], 2, max) >= 0
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

  # Get the connected components of the critical sets
  crits      = which(crit.set$minus)
  cut_points = which(diff(crits) != 1)
  cc.minus   = list()
  if(!identical(cut_points, integer(0))){

  }else{
    cc.minus[[1]] = crits
  }

  list(subI = subI, inter = subI.cit)
}

#' This function computes from a boolean vector the indices of the connected
#' components
#'
#' @param mask String name of the implemented method to
#' estimate the quantile function
#' @return A list containing for each entry the indices of the i-th
#'         connected components
#' @export
get_cc <- function(mask){
  ind_ccs    = which(mask)
  cut_points = which(diff(ind_ccs) != 1)
  ccs   = list()
  if(!identical(cut_points, integer(0))){
    cut_points = c(0, cut_points, length(ind_ccs))
    for(l in 1:(length(cut_points)-1)){
      ccs[[l]] = ind_ccs[(cut_points[l]+1):cut_points[l+1]]
    }
  }else{
    ccs[[1]] = ind_ccs
  }
  return(ccs)
}


#' This function computes from a boolean vector the indices of the connected
#' components
#'
#' @param mask String name of the implemented method to
#' estimate the quantile function
#' @return A list containing for each entry the indices of the i-th
#'         connected components
#' @export
knots2indices <- function(knots, x){
  inds = list()
  for(l in 1:(length(knots)-1)){
    inds[[l]] <- which(knots[l] <= x & x <= knots[l+1])
  }
  return(inds)
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
  arguments <- c(as.list(environment()), list(...))
  Narguments <- names(arguments)

  # Initialize the q_method list, which will be the output
  q.method <- list(
    name = arguments$name
  )

  # Fill the q_method list for different quantile function estimation
  # methods
  if(q.method$name %in% c("mboot", "fair.mboot")){
    # Use default values for the multiplier bootstrap, if
    # the values are not provided in the input arguments.
    if(!("Boottype" %in% Narguments)){
      q.method$Boottype = "t"
    }else{
      q.method$Boottype = arguments$Boottype
    }
    if(!("weights" %in% Narguments)){
      weights = "rademacher"
    }else{
      q.method$weights = arguments$weights
    }
    if(!("Mboots" %in% Narguments)){
      q.method$Mboots = 5e3
    }else{
      q.method$Mboots = arguments$Mboots
    }
    if(!("R" %in% Narguments)){
      stop("You need to provide residuals R.")
    }else{
      q.method$Mboots = arguments$R
    }

    # Add the fairness parameters if the bootstrap
    # should be fair
    if(q.method$name == "fair.mboot"){
      q.method$fair.intervals = arguments$fair.intervals
      if(!("type" %in% Narguments)){
        q.method$fair.type = "linear"
      }else{
        q.method$fair.type = arguments$fair.type
      }

      if(!("maxIter" %in% Narguments)){
        q.method$fair.niter = 10
      }else{
        q.method$fair.niter = arguments$fair.niter
      }
    }
  }else if(q.method$name == "gKR_t"){
    if(!("df" %in% Narguments)){
      q.method$df = NULL
    }else{
      q.method$df = df
    }
    if(!("knots" %in% Narguments)){
      q.method$knots = c(0,1)
    }else{
      q.method$knots = knots
    }
    q.method$Nknots = length(q.method$knots) - 1
    if(!("tau.est" %in% Narguments)){
      q.method$tau.est = tau_est
    }else{
      q.method$tau.est = tau.est
    }
    if(!("I_weights" %in% Narguments)){
      q.method$I_weights = rep(1/(length(q.method$knots) - 1), length(q.method$knots) - 1)
    }else{
      q.method$I_weights = I_weights
    }
    if(!("alpha_up" %in% Narguments)){
      q.method$alpha_up = alpha*(length(q.method$knots)-1)
    }else{
      q.method$alpha_up = alpha_up
    }
    if(!("maxIter" %in% Narguments)){
      q.method$maxIter = 0
    }else{
      q.method$maxIter = maxIter
    }
  }

  if("print.coverage" %in% Narguments){
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
  arguments <- c(as.list(environment()), list(...))
  Narguments <- names(arguments)

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

#' This function generates a list with the required input for different
#' methods to estimate the quantile function q of a SCoPE set.
#'
#' @param name String name of the implemented method to
#' estimate the quantile function
#' @param ... Optional parameters for the quantile estimation method
#' @return Scale field
#' @export
adapt_knotes_weights <- function(x, crit.set, knots, I_weights, vol = NULL){
  # x = seq(0,1,length.out=100)
  # dx = x[2]
  # crit.set <- list()
  # crit.set$minus = c(rep(T, 10), rep(F, 70), rep(T, 10), rep(F, 10))
  # crit.set$plus = c(rep(F, 5), rep(T, 3), rep(F, 56), T,rep(F, 4), rep(T, 8), rep(F, 3), rep(T, 2), rep(F, 3), rep(T, 13), rep(F, 2) )
  # knots = seq(0,1,length.out=5)
  # I_weights = rep(1/(length(knots) - 1), length(knots) - 1)
  #
  # plot(x, crit.set$minus,  col = 4, pch = 3)
  # points(x, crit.set$plus, col = 2, pch = 4)
  # abline(v=knots)
  #


  dx = diff(x)[1]

  s     = sub.intervals(x, knots, crit.set)
  subI  = s$subI
  inter = s$inter

  if(!all(inter)){
    I_weights = I_weights + sum(I_weights[!inter]) / sum(inter)
    I_weights[!inter] = 0
  }
  # Get the connected components of the union of the critical set
  ccs = get_cc(crit.set$minus | crit.set$plus)

  # Find additional knots
  subknots <- list()
  for(k in 1:length(I_weights)){
    tt = lapply(ccs, function(l) intersect(subI[[k]], l))

    ints_k <- NULL
    for(l in 1:length(tt)){
      if(length(tt[[l]]) != 0){
        ints_k <- c(ints_k, l)
      }
    }

    if(length(ints_k) >= 1){
      subknots[[k]] <- matrix(NaN, 3, length(ints_k))
      # Get the indices of the boundary of the subinterval
      inds_bdry_k = range(subI[[k]])
      # Go through the connected components within the subinterval
      for(ll in 1:length(ints_k)){
        # Find the knots which needs to be added and compute the volume
        # of the component within the subinterval
        if(length(tt[[ints_k[ll]]]) > 1){
          # Get the index range and the correspondig grid values
          inds_cc_ll = range(tt[[ints_k[ll]]])
          add_knots = x[inds_cc_ll]
          # remove the index
          if(inds_cc_ll[1] == inds_bdry_k[1] &&
             inds_cc_ll[2] == inds_bdry_k[2]){
            add_knots[1] = knots[k]
            add_knots[2] = knots[k+1]
            vol = diff(add_knots)
          }else if(inds_cc_ll[1] == inds_bdry_k[1]){
            add_knots[1] = knots[k]
            add_knots[2] = add_knots[2] + dx/2
            vol = diff(add_knots)
          }else if(inds_cc_ll[2] == inds_bdry_k[2]){
            add_knots[2] = knots[k+1]
            add_knots[1] = add_knots[1] - dx/2
            vol = diff(add_knots) + dx/2
          }else{
            add_knots = add_knots + c(-1, 1) * dx/2
            # Compute the volume
            vol = diff(add_knots)
          }
        }else{
          # treat single points as if they would be of size dx
          if(tt[[ints_k[ll]]] == 1){
            add_knots = c(1,1)*x[tt[[ints_k[ll]]]] + dx/2
            vol = dx/2
          }else if(tt[[ints_k[ll]]] == length(x)){
            add_knots = c(1,1)*x[tt[[ints_k[ll]]]] - dx/2
            vol = dx/2
          }else{
            add_knots = x[tt[[ints_k[ll]]]] + c(-1, 1) * dx/2
            vol = dx
          }
        }
        # save the knots and the volume for the ll'th cc
        subknots[[k]][,ll] <- c(add_knots, vol)
      }
      # Get the relative volume
      subknots[[k]][3,] <- subknots[[k]][3,] / sum(subknots[[k]][3,])
    }else{
      subknots[[k]] <- NA
    }
  }

  # Get the new knots and the new weights
  I_weights_new <- knots_new <- c()
  for(k in 1:length(I_weights)){
    if(!all(is.na(subknots[[k]]))){
      if(subknots[[k]][1,1] != knots[k] ){
        I_weights_new <- c(I_weights_new, 0)
        knots_new <- c(knots_new, knots[k])
      }
      for(i in 1:dim(subknots[[k]])[2] ){
        I_weights_new <- c(I_weights_new, c(I_weights[k]*subknots[[k]][3,i], 0) )
        knots_new     <- c(knots_new, subknots[[k]][1:2,i])
      }
      if(subknots[[k]][2, dim(subknots[[k]])[2]] == knots[k+1]){
        I_weights_new <- I_weights_new[-length(I_weights_new)]
        knots_new     <- knots_new[-length(knots_new)]
      }
    }else{
      I_weights_new <- c(I_weights_new, I_weights[k])
      knots_new <- c(knots_new, knots[k])
    }
  }
  knots_new <- c(knots_new, knots[length(knots)])

  # Iweights_fun = approxfun(x = new$knots[-length(new$knots)], y = new$I_weights,
  #                          method = "constant", rule = 2)
  # abline(v = new$knots, col = 3)
  # lines(x, Iweights_fun(x), col=2)

  # #-----------------------------------------------------------------------------
  # # Start half the weights if crit sets intersect
  # half    = (crit.set$minus == 1) & (crit.set$plus == 1)
  # if(any(half)){
  #   cc_half = get_cc(half)
  #
  #   knots     <- knots_new
  #   I_weights <- I_weights_new
  #
  #   for(l in 1:length(cc_half)){
  #     # Get the indices of the boundary of
  #     inds = range(cc_half[[l]])
  #     if(any(abs(x[inds[1]] - knots) <= dx/2 * (1 + 1e-8)) &
  #        any(abs(x[inds[2]] - knots) <= dx/2 * (1 + 1e-8))){
  #       # If both endpoints of the component agrees with
  #       # knots from the partition than half all of them between
  #       # these values
  #       ia = which(abs(x[inds[1]] - knots) <= dx/2 * (1 + 1e-8))[1]
  #       ie = which(abs(x[inds[2]] - knots) <= dx/2 * (1 + 1e-8))[1]
  #       for(i in ia:ie){
  #         I_weights[i] <- I_weights[i] / 2
  #       }
  #
  #     }else if(any(abs(x[inds[1]] - knots) <= dx/2 * (1 + 1e-8))){
  #       i = which(abs(x[inds[1]] - knots) <= dx/2 * (1 + 1e-8))
  #
  #       while( knots[i] < x[inds[2]] ){
  #         I_weights[i] <- I_weights[i] / 2
  #         i = i+1
  #       }
  #
  #       # Add a new knot and give the interval afterwards the same weight
  #       # as before but not modificated by 2
  #       knots     <- append(knots, x[inds[2]] + dx/2, after = i-1)
  #       I_weights <- append(I_weights, I_weights[i-1] * 2, after = i-1)
  #
  #     }else if(any(abs(x[inds[2]] - knots) <= dx/2 * (1 + 1e-8))){
  #       i = which(x[inds[1]] - knots < 0)[1] - 1
  #       # Add a new knot and give the interval afterwards the same weight
  #       # as before but not modificated by 2
  #       knots     <- append(knots, x[inds[1]] - dx/2, after = i)
  #       I_weights <- append(I_weights, I_weights[i] / 2, after = i)
  #
  #       i = i + 2
  #       while( knots[i] < x[inds[2]] ){
  #         I_weights[i] <- I_weights[i] / 2
  #         i = i+1
  #       }
  #
  #     }else{
  #       i = which(x[inds[1]] - knots < 0)[1] - 1
  #       # Add a new knot and give the interval afterwards the same weight
  #       # as before but not modificated by 2
  #       knots     <- append(knots, x[inds[1]] - dx/2, after = i)
  #       I_weights <- append(I_weights, I_weights[i] / 2, after = i)
  #
  #       i = i + 2
  #       while( knots[i] < x[inds[2]] ){
  #         I_weights[i] <- I_weights[i] / 2
  #         i = i+1
  #       }
  #
  #       # Add a new knot and give the interval afterwards the same weight
  #       # as before but not modificated by 2
  #       knots     <- append(knots, x[inds[2]] + dx/2, after = i-1)
  #       I_weights <- append(I_weights, I_weights[i-1] * 2, after = i-1)
  #     }
  #   }
  # }

  # plot(x, crit.set$minus,  col = 4, pch = 3)
  # points(x, crit.set$plus, col = 2, pch = 4)
  # abline(v = knots, col = 3)
  #
  # Iweights_fun = approxfun(x = knots[-length(knots)], y = I_weights,
  #                          method = "constant", rule = 2)
  # lines(x, Iweights_fun(x), col=2)
  #

  return(list(knots = knots, I_weights = I_weights))
}
