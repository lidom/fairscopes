#------------------------------------------------------------------------------#
#                                                                              #
#     Predefined functions to estimate to estimate the tau parameter           #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' Approximates the true tau parameter from a given covariance
#' function using numeric derivatives.
#'
#' @param x (numeric vector) a vector for evaluation
#' @param cov_fun (function f(x,y,...)) a covariance function.
#' @param h (numeric) the stepsize for numerical approximation of the derivative
#' @return
#' A vector containing the numerical approximation of tau(x).
#' @export
true_tau <- function( x, cov_fun, dx = 1e-8, ...){
  x_p = x + dx
  x_p[length(x)] = x_p[length(x)-1]
  x_m = x - dx
  x_m[1] = x_m[2]

  return(sqrt((vapply(1:length(x), function(l) cov_fun(x_p[l], x_p[l], ...), FUN.VALUE = 0.1) -
      2*vapply(1:length(x), function(l) cov_fun(x_p[l], x_p[l], ...), FUN.VALUE = 0.1) +
        vapply(1:length(x), function(l) cov_fun(x_m[l], x_m[l], ...), FUN.VALUE = 0.1) ) / 4 / dx^2))
}


#' Estimates tau parameter by interpolation of the data by
#' natural splines and then numerical differentiation.
#' There is the option to smooth tau afterwards using a
#' smoothing spline.
#'
#' @param x a vector for evaluation
#' @param cov_fun a covariance function.
#' @param h stepsize for approximation of the derivative
#' @param df Default NULL. If not null, the estimated tau
#' function from interpolation is smoothed using a smoothing
#' spline. Here df is explained in smooth.spline().
#' @param ... Options of smooth.spline()
#' @return
#' The approximation of tau at x.
#' @export
tau_est <- function(R, x, df = NULL, ...){
  dR <- apply( R, 2,
               FUN = function( yy ){
                 # Interpolate the function
                 fn <- stats::splinefun( x = x,
                                         y = yy,
                                         method = "natural" )
                 # Obtain the derivative of the function
                 pracma::fderiv( f = fn,
                                 x = x,
                                 n = 1,
                                 method = "central")
               } )

  sd.dR = apply( dR , 1, stats::sd )
  if(is.null(df)){
    return(approxfun(x, y = sd.dR, method = "linear", yleft = 0, yright = 0))
  }else{
    test = smooth.spline(x, y = sd.dR, df = df, ...)
    return(approxfun(x, y = test$y, method = "linear", yleft = 0, yright = 0))
  }
}


# #' Estimates tau parameter by interpolation of the data by
# #' natural splines and then numerical differentiation.
# #' There is the option to smooth tau afterwards using a
# #' smoothing spline.
# #'
# #' @param x a vector for evaluation
# #' @param cov_fun a covariance function.
# #' @param h stepsize for approximation of the derivative
# #' @param df Default NULL. If not null, the estimated tau
# #' function from interpolation is smoothed using a smoothing
# #' spline. Here df is explained in smooth.spline().
# #' @param ... Options of smooth.spline()
# #' @return
# #' The approximation of tau at x.
# #' @export
# tau_est_basis <- function(R, x, nbasis = 40, df = NULL, ...){
#   #
#   N = dim(R)[2]
#
#   #
#   Basis = fda::create.bspline.basis(rangeval = range(x), nbasis = nbasis)
#
#   # Get fd object from Basis and its derivatives
#   Basis.fd  = fda::fd(diag(rep(1, nbasis)), Basis)
#   dBasis.fd = fda::deriv.fd(Basis.fd)
#
#   fd1 <- fda::smooth.basis(argvals = x, y =  R, Basis)$fd
#
#   # par(mfrow=c(1,2))
#   # matplot(x, R, type = "l")
#   # matplot(x, fda::eval.fd(x , fd1), type = "l")
#
#   # get the empirical estimate of the asymptotic covariance estimate in basis
#   C = (N-1)/N * cov(t(fd1$coefs))
#
#   # Estimate LKC L1
#   bx     = fda::eval.fd(seq(0, 1, length.out = Nx), dBasis.fd)
#   sqdvar = apply(bx, 1, function(r){ sqrt(t(r) %*% C %*% r)} )
#
#   if(is.null(df)){
#     return(approxfun(x, y = sd.dR, method = "linear", yleft = 0, yright = 0))
#   }else{
#     test = smooth.spline(x, y = sd.dR, df = df, ...)
#     return(approxfun(x, y = test$y, method = "linear", yleft = 0, yright = 0))
#   }
# }
