#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to compute SCoPE sets
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
# - fairSCB()
# - fairSCB_var()
# - plot_SCoPES()
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' This functions computes fair SCBs given an estimator an estimated
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
fairSCB <- function(alpha, hatmu, hatrho, tN,
                    x = seq(0, 1, length.out = length(hatmu)),
                    q.method, type = "two-sided", mu = NULL, subI = NULL){

  #---------------------------------------------------------------------------
  # Estimate the quantile funcion q.
  if(type == "two-sided"){
    crit.set = list(minus = rep(TRUE, length(hatmu)),
                    plus  = rep(TRUE, length(hatmu)))
  }else{
    crit.set = list(minus = rep(TRUE, length(hatmu)),
                    plus  = rep(FALSE, length(hatmu)))
  }


  fair.q = NULL
  if(q.method$name == "mboot"){
    q = MultiplierBootstrapSplit(alpha   = alpha,
                                 R       = q.method$R,
                                 minus   = crit.set$minus,
                                 plus    = crit.set$plus,
                                 Mboots  = q.method$Mboots,
                                 method  = q.method$Boottype,
                                 weights = q.method$weights)$q
    if(is.infinite(q)){
      q = 0
    }

    width           <- matrix(0, nrow = 1, ncol = 2)
    colnames(width) <- c("L1", "L2")
    rownames(width) <- c("u")
    width[1,]       <- c(q,q) * diff(range(x))

    u <- Vectorize(function(x) q)
    q <- u(x)

    # Preapare the output function
    out = list()

  }else if(q.method$name == "fair.mboot"){
    samples = MultiplierBootstrapSplit(alpha   = alpha,
                                       R       = q.method$R,
                                       minus   = crit.set$minus,
                                       plus    = crit.set$plus,
                                       Mboots  = q.method$Mboots,
                                       method  = q.method$type,
                                       weights = q.method$weights)$samples
    if( !(all(!crit.set$minus) && all(!crit.set$plus)) ){
      fair.q = fair_quantile_boot(alpha  = alpha,
                                  x = x,
                                  samples   = samples,
                                  crit.set  = crit.set,
                                  knots     = q.method$knots,
                                  I_weights = q.method$I_weights,
                                  type      = q.method$type,
                                  alpha_up  = q.method$alpha_up,
                                  maxIter   = q.method$maxIter)
      u = fair.q$u
      width           <- matrix(0, nrow = 1, ncol = 2)
      colnames(width) <- c("L1", "L2")
      rownames(width) <- c("u")

      width[1,] <- c(integrate(function(t) abs(u(t)),   lower = 0, upper = 1)$value,
                     sqrt(integrate(function(t) abs(u(t))^2, lower = 0, upper = 1)$value))

      q = u(x)
    }else{
      fair.q = NULL
      q = rep(0, length(x))
    }

    # Preapare the output function
    out = list()

  }else if(q.method$name == "KRF_fair"){
    if(type == "two-sided"){
      alpha = alpha / 2
    }

    u = fair_quantile_KRF(alpha     = alpha,
                          type      = q.method$type,
                          tau       = q.method$tau,
                          df        = q.method$df,
                          knots     = q.method$knots,
                          I_weights = q.method$I_weights,
                          u.type    = q.method$u.type,
                          alpha_up  = q.method$alpha_up,
                          maxIter   = q.method$maxIter,
                          tol       = alpha / 100)$u

    width           <- matrix(0, nrow = 1, ncol = 2)
    colnames(width) <- c("L1", "L2")
    rownames(width) <- c("u")

    width[1,] <- c(integrate(function(t) abs(u(t)),   lower = 0, upper = 1)$value,
                   sqrt(integrate(function(t) abs(u(t))^2, lower = 0, upper = 1)$value))

    q = u(x)

    # Preapare the output function
    out <- list()

  }else if(q.method$name == "KRF_nonfair"){
    if(type == "two-sided"){
      alpha = alpha / 2
    }

    q = RFT::GKFthreshold( alpha = alpha,
                           LKC = c(1, integrate_save(q.method$tau,
                                                xlims = c(x[1], x[length(x)]))),
                           type = q.method$type,
                           df   = q.method$df,
                           interval = c(0, 100) )$threshold

    width           <- matrix(0, nrow = 1, ncol = 2)
    colnames(width) <- c("L1", "L2")
    rownames(width) <- c("u")
    width[1,]       <- c(q,q) * diff(range(x))

    u <- Vectorize(function(x) q)
    q <- u(x)

    # Preapare the output function
    out <- list()

  }else if(q.method$name == "KRF_min"){
    if(type == "two-sided"){
      alpha = alpha / 2
    }

    K = length(q.method$knots)-1
    # Get the starting value
    u0 = RFT::GKFthreshold( alpha = alpha,
                           LKC = c(1, integrate_save(q.method$tau,
                                                     xlims = c(x[1], x[length(x)]))),
                           type = q.method$type,
                           df   = q.method$df,
                           interval = c(0, 100) )$threshold

    # Find the
    thresh <- quantile_KRF(x0 = c(u0, rep(0,K)),
                           alpha,
                           tau   = q.method$tau,
                           type  = q.method$type,
                           df    = q.method$df,
                           knots = q.method$knots,
                           alpha_weights = q.method$I_weights)

    q <- thresh$u(x)
    u <- thresh$u

    width <- thresh$width

    out <- list(optim = thresh$res)
    # out$constraints_check = thresh$constraints_check
  }

  SCB = data.frame(
    "x"   = x,
    "low" = hatmu - tN*hatrho*q,
    "est" = hatmu,
    "up"  = hatmu + tN*hatrho*q,
    "q"   = q
  )

  if(type == "low"){
    SCB = SCB[-4]
  }else if(type == "up"){
    SCB = SCB[, -2]
  }

  # Output
  if(is.null(mu)){
    return(c(out,
             list(SCB = SCB, width = width, u = u)))
  }else{
    loc.cov = abs(mu - hatmu) <= q*tN*hatrho
    SCB$mu  = mu
    # Note that the coverage only works for two-sided!
    return(c(out,
             list(SCB = SCB, width = width, u = u,
                loc.cov  = loc.cov,
                glob.cov = all(loc.cov))))
  }
}


#' This functions computes fair SCBs given an estimator an estimated
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
fairSCB_var <- function(alpha, hatvar,
                        x = seq(0, 1, length.out = length(hatvar)),
                        q.method, type = "two-sided", true_var = NULL,
                        subI = NULL){
  #---------------------------------------------------------------------------
  # Estimate the quantile funcion q.
  fair.q = NULL
  if(q.method$name == "KRF_fair"){
    if(type == "two-sided"){
      lb = fair_quantile_KRF(alpha = alpha/2,
                             type  = "chi2",
                             tau   = q.method$tau,
                             df        = q.method$df,
                             knots     = q.method$knots,
                             I_weights = q.method$I_weights,
                             u.type      = "linear",
                             crossings = "up",
                             lower.tail = FALSE,
                             alpha_up  = q.method$alpha/2,
                             maxIter   = q.method$maxIter,
                             tol = alpha/2 / 100)
      ub = fair_quantile_KRF(alpha = alpha/2,
                             type  = "chi2",
                             tau   = q.method$tau,
                             df         = q.method$df,
                             knots      = q.method$knots,
                             I_weights  = q.method$I_weights,
                             u.type     = "linear",
                             crossings  = "down",
                             lower.tail = TRUE,
                             alpha_up   = q.method$alpha/2,
                             maxIter    = q.method$maxIter,
                             tol = alpha/2 / 100)

      width           <- matrix(0, nrow = 1, ncol = 2)
      colnames(width) <- c("L1", "L2")
      rownames(width) <- c("u")

      width[1,] <- c(integrate(function(t) abs(1/ub$u(t) - 1/lb$u(t)),   lower = 0, upper = 1)$value,
                     sqrt(integrate(function(t) abs(1/ub$u(t) - 1/lb$u(t))^2, lower = 0, upper = 1)$value))

    }else if(type == "low"){
      lb = fair_quantile_KRF(alpha = alpha,
                             type  = "chi2",
                             tau   = q.method$tau,
                             df        = q.method$df,
                             knots     = q.method$knots,
                             I_weights = q.method$I_weights,
                             u.type      = "linear",
                             crossings   = "up",
                             lower.tail  = FALSE,
                             alpha_up  = q.method$alpha,
                             maxIter   = q.method$maxIter,
                             tol = alpha / 100)
      ub <- list()
      ub$ub <- Vectorize(function(t){Inf})

      SCB   = cbind(q.method$df*hatvar/lb$u(x), rep(Inf, length(x)))
      width = NULL
    }else if(type == "up"){
      ub = fair_quantile_KRF(alpha = alpha,
                             type  = "chi2",
                             tau   = q.method$tau,
                             df         = q.method$df,
                             knots      = q.method$knots,
                             I_weights  = q.method$I_weights,
                             u.type     = "linear",
                             crossings  = "low",
                             lower.tail = TRUE,
                             alpha_up   = q.method$alpha,
                             maxIter    = q.method$maxIter,
                             tol = alpha / 100)

      lb <- list()
      lb$lb <- Vectorize(function(t){0})

      width = NULL
    }
    out <- list()
  }else if(q.method$name == "KRF_nonfair"){
    if(type == "two-sided"){
      alpha = alpha / 2

      ub =  fair_quantile_KRF(alpha = alpha,
                               type  = "chi2",
                               tau   = q.method$tau,
                               df         = q.method$df,
                               knots      = range(x),
                               I_weights  = 1,
                               u.type     = "linear",
                               crossings  = "low",
                               lower.tail = TRUE,
                               alpha_up   = q.method$alpha,
                               maxIter    = 0,
                               tol = alpha / 100)

      lb = fair_quantile_KRF(alpha = alpha,
                              type  = "chi2",
                              tau   = q.method$tau,
                              df         = q.method$df,
                              knots      = range(x),
                              I_weights  = 1,
                              u.type     = "linear",
                              crossings  = "up",
                              lower.tail = FALSE,
                              alpha_up   = q.method$alpha,
                              maxIter    = 0,
                              tol = alpha / 100)

      width           <- matrix(0, nrow = 1, ncol = 2)
      colnames(width) <- c("L1", "L2")
      rownames(width) <- c("u")
      width[1,]       <- c(1/ub$u(0)-1/lb$u(0),1/ub$u(0)-1/lb$u(0)) * diff(range(x))
    }
    # Preapare the output function
    out <- list()

  }else if(q.method$name == "KRF_min"){
    if(type == "two-sided"){
      alpha = alpha / 2
    }

    K = length(q.method$knots)-1

    # Get the starting value
    ub0 =  fair_quantile_KRF(alpha = alpha,
                            type  = "chi2",
                            tau   = q.method$tau,
                            df         = q.method$df,
                            knots      = range(x),
                            I_weights  = 1,
                            u.type     = "linear",
                            crossings  = "low",
                            lower.tail = TRUE,
                            alpha_up   = q.method$alpha,
                            maxIter    = 0,
                            tol = alpha / 100)

    lb0 = fair_quantile_KRF(alpha = alpha,
                           type  = "chi2",
                           tau   = q.method$tau,
                           df         = q.method$df,
                           knots      = range(x),
                           I_weights  = 1,
                           u.type     = "linear",
                           crossings  = "up",
                           lower.tail = FALSE,
                           alpha_up   = q.method$alpha,
                           maxIter    = 0,
                           tol = alpha / 100)

    # Find the
    thresh_ub <- quantile_KRF(x0 = c(ub0$u(0), rep(0,K)),
                           alpha = alpha,
                           tau   = q.method$tau,
                           type  = q.method$type,
                           df    = q.method$df,
                           knots = q.method$knots,
                           alpha_weights = q.method$I_weights,
                           crossings  = "down",
                           lower.tail = TRUE,
                           MAX = TRUE)
    thresh_lb <- quantile_KRF(x0 = c(lb0$u(0), rep(0,K)),
                              alpha,
                              tau   = q.method$tau,
                              type  = q.method$type,
                              df    = q.method$df,
                              knots = q.method$knots,
                              alpha_weights = q.method$I_weights,
                              crossings  = "up",
                              lower.tail = FALSE,
                              MAX = FALSE)

    width           <- matrix(0, nrow = 1, ncol = 2)
    colnames(width) <- c("L1", "L2")
    rownames(width) <- c("u")

    ub <- thresh_ub
    lb <- thresh_lb

    width[1,] <- c(integrate(function(t) abs(1/ub$u(t) - 1/lb$u(t)),   lower = 0, upper = 1)$value,
                   sqrt(integrate(function(t) abs(1/ub$u(t) - 1/lb$u(t))^2, lower = 0, upper = 1)$value))

    out <- list(optim = list(lb = thresh_lb$res, ub = thresh_ub$res))
    # out$constraints_check = thresh$constraints_check
  }

  SCB = data.frame(
    "x"   = x,
    "low" = q.method$df * hatvar / lb$u(x),
    "est" =  hatvar,
    "up"  = q.method$df * hatvar / ub$u(x),
    "qlow" = lb$u(x),
    "qup"  = ub$u(x)
  )

  if(type == "low"){
    SCB = SCB[-4]
  }else if(type == "up"){
    SCB = SCB[, -2]
  }

  # Output
  if(is.null(true_var)){
      return(c(out,
             list(SCB = SCB, width = width, ub = ub, lb = lb)))
  }else{
      loc.cov = (SCB[, "low"] <= true_var) & (true_var <= SCB[, "up"])

      SCB$var  = true_var
      # Note that the coverage only works for two-sided!
      return(c(out,
               list(SCB = SCB, width = width, ub = ub, lb = lb,
                    loc.cov  = loc.cov,
                    glob.cov = all(loc.cov))))
  }
}


#' This functions computes fair SCBs given an estimator an estimated
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
fairInferencePtwLinModel <-
  function(x, Y, X, W = NULL, sigma = NULL,
           alpha,  q.method, type = "two-sided", mu = NULL, subI = NULL){

    #---------------------------------------------------------------------------
    # Estimate the quantile funcion q.
    if(type == "two-sided"){
      crit.set = list(minus = rep(TRUE, length(hatmu)),
                      plus  = rep(TRUE, length(hatmu)))
    }else{
      crit.set = list(minus = rep(TRUE, length(hatmu)),
                      plus  = rep(FALSE, length(hatmu)))
    }


    fair.q = NULL
    if(q.method$name == "mboot"){
      q = MultiplierBootstrapSplit(alpha   = alpha,
                                   R       = q.method$R,
                                   minus   = crit.set$minus,
                                   plus    = crit.set$plus,
                                   Mboots  = q.method$Mboots,
                                   method  = q.method$Boottype,
                                   weights = q.method$weights)$q
      if(is.infinite(q)){
        q = 0
      }
    }else if(q.method$name == "fair.mboot"){
      samples = MultiplierBootstrapSplit(alpha   = alpha,
                                         R       = q.method$R,
                                         minus   = crit.set$minus,
                                         plus    = crit.set$plus,
                                         Mboots  = q.method$Mboots,
                                         method  = q.method$Boottype,
                                         weights = q.method$weights)$samples
      if( !(all(!crit.set$minus) && all(!crit.set$plus)) ){
        fair.q = fair_quantile_boot(alpha  = alpha,
                                    x = x,
                                    samples = samples,
                                    crit.set  = crit.set,
                                    knots     = q.method$knots,
                                    I_weights = q.method$I_weights,
                                    type      = q.method$type,
                                    alpha_up  = q.method$alpha_up,
                                    maxIter   = q.method$maxIter)
        q = fair.q$u(x)
      }else{
        fair.q = NULL
        q = rep(0, length(x))
      }

    }else if(q.method$name == "KR_t"){
      if(type == "two-sided"){
        alpha = alpha / 2
      }

      q = fair_quantile_KRF(alpha = alpha,
                            type  = "t",
                            tau   = q.method$tau,
                            df        = q.method$df,
                            knots     = q.method$knots,
                            I_weights = q.method$I_weights,
                            u.type      = "linear",
                            alpha_up  = q.method$alpha_up,
                            maxIter   = q.method$maxIter,
                            tol = alpha / 100)$u
      q = q(x)
    }else if(q.method$name == "KR_t_nonfair"){
      if(type == "two-sided"){
        alpha = alpha / 2
      }

      q = RFT::GKFthreshold( alpha = alpha,
                             LKC = c(1, integrate(q.method$tau,
                                                  lower = x[1],
                                                  upper = x[length(x)])$val),
                             type = "t",
                             df = q.method$df,
                             interval = c(0, 100) )$threshold
    }else if(q.method$name == "KR_z_nonfair"){
      if(type == "two-sided"){
        alpha = alpha / 2
      }

      q = RFT::GKFthreshold( alpha = alpha,
                             LKC = c(1, integrate(q.method$tau,
                                                  lower = x[1],
                                                  upper = x[length(x)])$val),
                             type = "z",
                             df = 0,
                             interval = c(0, 100) )$threshold
    }

    SCB = cbind(hatmu - tN*hatrho*q, hatmu + q*tN*hatrho)
    rownames(SCB) <- x
    colnames(SCB) <- c("low", "up")

    if(type == "low"){
      SCB = SCB[, "low"]
    }else if(type == "low"){
      SCB = SCB[, "up"]
    }

    # Output
    if(is.null(mu)){
      return(list(SCB = SCB, q = q))
    }else{
      loc.cov = abs(mu - hatmu) <= q*tN*hatrho

      # Note that the coverage only works for two-sided!
      return(list(SCB = SCB,
                  q = q,
                  loc.cov  = loc.cov,
                  glob.cov = all(loc.cov)))
    }
  }
