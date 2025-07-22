#------------------------------------------------------------------------------#
#                                                                              #
#     Constraints and loss functions                                           #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#
# Contained functions:
#      - L1L2_loss() (removed)
#      - KRF_constraint
#------------------------------------------------------------------------------#
# Developer notes:
# The losses depend on a basis function expansion. At the moment the package
# fda is used for it.
#------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------
# Loss functions
#-------------------------------------------------------------------------------
#' Computes the L1 and L2 norm of a function given in a basis expansion
#'
#' @param coef array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param basis an basis object from fda package
#' @return vector of length 2 containing the L1 norm and square of the L2 norm
#' of the function defined by coef and basis
#' @export
L1L2_loss <- function(coef, basis) {
  # initialize the output function
  out <- c(0,0)
  names(out) <- c("L1", "L2")

  # get the functional data object
  u <- fd(coef = coef, basisobj = basis)

  # compute L1 and L2 norm of u
  out[1] = integrate(function(t) abs(eval.fd(t, u)),
                     lower = 0, upper = 1)$value
  out[2] = sqrt(inprod(u, u))

  # Output the two norms
  return(out)
}


#-------------------------------------------------------------------------------
# Constraints
#-------------------------------------------------------------------------------
#' Computes the L1 and L2 norm of a function given in a basis expansion
#'
#' @param coef array of dimension K x N containing N-realizations of
#'  a random field over a 1-dimensional domain.
#' @param basis an basis object from fda package
#' @return vector of length 2 containing the L1 norm and square of the L2 norm
#' of the function defined by coef and basis
#' @export
quantile_KRF <- function(alpha, tau, basis, type = "t", df = NULL,
                         loss_type     = c(0, 1, 0.2, 0),
                         knots         = c(0, 1),
                         alpha_weights = rep(1 / (length(knots) - 1), length(knots) - 1),
                         print = TRUE){
  # Initialize output
  summary <- list()
  # Get the necessary parameters from the input
  #-----------------------------------------------------------------------------
  # Get the amount of intervals
  NI = length(knots)-1
  # Get the correct KRF formula
  if(type == "t"){
    KRF = KacRice_t
  }

  # Get the starting value for the minimization
  #-----------------------------------------------------------------------------
  KRFu_up = function(y){
    KRF(df = df,
        tau = tau,
        u   = Vectorize(function(x) y),
        du  = Vectorize(function(x) 0),
        x   = knots[c(1, NI+1)],
        crossings = "up", t0 = 1, sigma = 1,
        lower.tail = FALSE, EC = TRUE) - alpha}
  tt = uniroot(f = KRFu_up, interval = c(0,10))

  # Currently assumes basis to be B-splines
  x0 <- c(1, rep(tt$root, nbasis))
  q0 <- fd(coef = x0[-1], basisobj = basis)

  # Define the width loss
  #-----------------------------------------------------------------------------
  WidthLoss <- function(x){
    # get the function u and its derivatives
    u   <- fd(coef = x[-1], basisobj = basis)
    du  <- deriv.fd(u, Lfd = 1)
    d2u <- deriv.fd(u, Lfd = 2)

    # initialize the output
    out = 0

    if(loss_type[1] != 0){
      # Funktion, die an t den absoluten Funktionswert zurückgibt
      f_abs <- function(t) abs(eval.fd(t, u))

      out = out + integrate(f_abs, lower = 0, upper = 1)$value
    }

    if(loss_type[2] != 0){
      out = out + inprod(u, u)
    }

    if(loss_type[3] != 0){
      out = out + loss_type[3] * inprod(du, du)
    }

    if(loss_type[4] != 0){
      out = out + loss_type[4] * inprod(d2u, d2u)
    }
    out
  }

  # Define the different constraints
  #-----------------------------------------------------------------------------
  # equality constraints
  KRF_constraint <- function(x) {
    q  <- fd(coef = x[-1], basisobj = basis)
    dq <- deriv.fd(q, Lfd = 1)
    u  = function(y){ as.vector(eval.fd(evalarg = y, fdobj = q)) }
    du = function(y){ as.vector(eval.fd(evalarg = y, fdobj = dq)) }

    if(NI == 1){
      constraints <- KRF(df, tau, u, du, x = c(knots[1], knots[NI+1]), t0 = 1) - alpha
    }else{
      # inner knots
      xx = knots[-c(1, NI+1)]

      constraints <- c(vapply(1:NI, function(i)
        KRF(df, tau, u, du, x = c(knots[i], knots[i+1]), t0 = 1) - x[1]*alpha_weights[i]*alpha,
        FUN.VALUE = pi),
        sum(pnorm(u(xx), lower.tail = FALSE)) - (x[1]-1) * alpha)
    }

    # Output the evalated constraints
    constraints
  }

  # Inequality constraint
  # Computes the inequality constraint on the parameter c, i.e. x[1].
  ineq_constraint <- function(x) {
    1 - x[1]
  }

  if( !is.null(loss_type) ){
    # Optimization
    #-----------------------------------------------------------------------------
    res <- nloptr(
      x0 = x0,
      eval_f = WidthLoss,
      eval_g_eq = KRF_constraint,
      eval_g_ineq = ineq_constraint,
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

    # Summarize results
    #-----------------------------------------------------------------------------
    summary[["results"]] <- res

    u_optim = res$solution
    q = fd(coef = u_optim[-1], basisobj = basis)
    dq = deriv.fd(q, Lfd = 1)
    u  = function(y){ as.vector(eval.fd(evalarg = y, fdobj = q)) }
    du = function(y){ as.vector(eval.fd(evalarg = y, fdobj = dq)) }

    summary[["u"]]  <- u
    summary[["du"]] <- du

    # Get the constraints for the optimized function to check whether they are satisfied
    if(NI > 1){
      constraints_check = round(c(KRF_constraint(u_optim),
                                  KRF(df, tau, u, du, x = knots[c(1, NI+1)], t0 = 1)), 4)
      names(constraints_check) <- c(paste("err:", 1:NI, sep=""), "ErrTot", "KRF")
    }else{
      constraints_check = round(c(KRF_constraint(u_optim),
                                  KRF(df, tau, u, du, x = knots[c(1, NI+1)], t0 = 1)), 4)
      names(constraints_check) <- c("ErrTot", "KRF")
    }
    summary[["constraints_check"]] <- constraints_check

    # Get the L1 and L2 indicators before and after minimization
    WidthSCB <- matrix(0, nrow = 3, ncol = 2)
    colnames(WidthSCB) <- c("L1", "L2")
    rownames(WidthSCB) <- c("start", "opt", "rel. improv.")
    WidthSCB[1,] <- L1L2_loss(x0[-1], basis = basis)
    WidthSCB[2,] <- L1L2_loss(res$solution[-1], basis = basis)
    WidthSCB[3,] <- round(100*(diff(WidthSCB)[1,]/WidthSCB[1,]),2)

    summary[["WidthSCB"]] <- WidthSCB

    # Plot the tau function and the optimized quantile
    if(print == TRUE){
      eval_pts = seq(knots[1], knots[NI+1], length.out=100)
      par(mfrow = c(1,2))
      plot(eval_pts, tau(eval_pts),
           main = "tau function", type = "l",
           xlab = "", ylab = "", lwd = 2)
      abline(v = knots, lty = 2, col = "grey70")

      plot(q, main = paste("Threshold Function, rel. gain: L1=", summary$WidthSCB[3,1],
                           ", L2=", summary$WidthSCB[3,2]), lwd = 2)
      abline(v = knots, lty = 2, col = "grey70")
      abline(h = tt$root, col = "red", lwd = 2)
      legend("topleft", legend = c("Start", "Optim."),
             col = c("red", "black"), lty = 1, lwd = 2,
             bg = "transparent",  bty = "n")
    }
  }else{
    u  <- Vectorize(function(x) tt$root)
    du <- Vectorize(function(x) 0)
    summary[["u"]]  <- u
    summary[["du"]] <- du

    # Get the constraints for the optimized function to check whether they are satisfied
    if(NI > 1){
      constraints_check = round(c(KRF_constraint(x0),
                                  KRF(df, tau, u, du, x = knots[c(1, NI+1)], t0 = 1)), 4)
      names(constraints_check) <- c(paste("err:", 1:NI, sep=""), "ErrTot", "KRF")
    }else{
      constraints_check = round(c(KRF_constraint(x0),
                                  KRF(df, tau, u, du, x = knots[c(1, NI+1)], t0 = 1)), 4)
      names(constraints_check) <- c("ErrTot", "KRF")
    }
    summary[["constraints_check"]] <- constraints_check

    WidthSCB <- matrix(0, nrow = 1, ncol = 2)
    colnames(WidthSCB) <- c("L1", "L2")
    rownames(WidthSCB) <- c("u")
    WidthSCB[1,] <- L1L2_loss(x0[-1], basis = basis)
    summary[["WidthSCB"]] <- WidthSCB
  }

  summary
}

