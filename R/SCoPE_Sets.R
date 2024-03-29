#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to compute SCoPE sets
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
# - hatLc()
# - hatUc()
# - SCoPES()
# - plot_SCoPES()
# - fairSCB()
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' This functions computes an lower confidence excursion set
#'
#' @inheritParams PreimageC_thickening
#' @param c (numeric vector) vector with the same
#' length as 'hatmu' containing the thresholding function
#' @param q (numeric vector) quantile vector for the SCoPE set has the same
#' length as 'hatmu'
#' @param mu (numeric vector) vector containing the true values of the target
#' function mu.
#' @param inclusion (string) if "inequal" the lower excursion set is open, if "equal"
#' the lowerer excursion set is closed.
#' @return  if is.null(mu) a boolean vector indicating which positions are
#' in the lower confidence excursion set. Otherwise a list with elements
#'  \itemize{
#'   \item hatUc boolean vector indicating which positions are
#' in the lower confidence excursion set
#'   \item contain boolean indicating whether the lower confidence set is contained in
#'   the true lower excursion set
#' }
#' @export
hatLc <- function(c, hatmu, hatsigma, tN, q, mu = NULL, inclusion = "inequal" ){
  if(inclusion == "inequal"){
    hatL = hatmu < c - q * tN * hatsigma
    if(!is.null(mu)){
      Lc = mu < c
      contain  = (hatL | Lc) == Lc
      detect  <- NA * (1:length(hatmu))
      detect[hatL] <- contain[hatL]
      return(list(hatLc = hatL, contain = contain, detect = detect))
    }else{
      return(hatL)
    }
  }else{
    hatU = hatUc(c = c, hatmu = hatmu, hatsigma = hatsigma,
                 tN = tN, q = q, mu = mu, inclusion = "inequal" )
    if(!is.null(mu)){
      return(list(hatLc = !hatU$hatUc, contain = hatU$contain))
    }else{
      return(!hatU)
    }
  }
}


#' This functions computes an upper confidence excursion set
#'
#' @inheritParams PreimageC_thickening
#' @inheritParams hatLc
#' @return  if is.null(mu) a boolean vector indicating which positions are
#' in the upper confidence excursion set. Otherwise a list with elements
#'  \itemize{
#'   \item hatUc boolean vector indicating which positions are
#' in the upper confidence excursion set
#'   \item contain boolean indicating whether the upper confidence set is contained in
#'   the true upper excursion set
#' }
#' @export
hatUc <- function(c, hatmu, hatsigma, tN, q, mu = NULL, inclusion = "inequal" ){
  if(inclusion == "inequal"){
    hatU = hatmu > c + q * tN * hatsigma
    if(!is.null(mu)){
      Uc = mu > c
      contain <- (hatU | Uc) == Uc
      detect  <- NA * (1:length(hatmu))
      detect[hatU] <- contain[hatU]
      return(list(hatUc = hatU, contain = contain, detect = detect))
    }else{
      return(hatU)
    }
  }else{
    hatL = hatLc(c = c, hatmu = hatmu, hatsigma = hatsigma,
                 tN = tN, q = q, mu = mu, inclusion = "inequal" )
    if(!is.null(mu)){
      return(list(hatUc = !hatL$hatLc, contain = hatL$contain))
    }else{
      return(!hatL)
    }
  }
}


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams PreimageC_thickening
#' @inheritParams hatLc
#' @param alpha (numeric) between 0 and 1. It will produce (1-alpha) SCoPES.
#' @param x (numeric vector) a vector of length hatmu which contains the
#' locations/times of the function mu etc, i.e., it enumerates discrete
#' samples from S. Default seq(0, 1, length.out=length(hatmu)).
#' @param q.method (list) ...
#' @param Preimage.method (list) ...
#' @param inclusion list with elements
#'  \itemize{
#'   \item L string either "inequal" or "equal" determining whether the lower
#'   excursion set uses < c or <=c.
#'   \item U string either "inequal" or "equal" determining whether the upper
#'   excursion set uses > c or >=c.
#' }
#' @return Standard error under the assumption the data is Gaussian
#' @export
SCoPES <- function(alpha, C, x = seq(0, 1, length.out = length(hatmu)),
                   hatmu, hatsigma, tN,
                   q.method,
                   Preimage.method,
                   inclusion = list(L = "inequal", U = "inequal"),
                   mu = NULL){
  #-----------------------------------------------------------------------------
  # Get standard values and catch wrong inputs

  # Get the number of levels and make C a column matrix if
  # it is just one function.
  num.levels = list()
  if(!is.null(C$minus)){
    if(is.null(dim(C$minus))){
      num.levels$minus = 1
      C$minus = t(t(C$minus))
    }else{
      num.levels$minus = dim(C$minus)[2]
    }
  }
  if(!is.null(C$plus)){
    if(is.null(dim(C$plus))){
      num.levels$plus = 1
      C$plus = t(t(C$plus))
    }else{
      num.levels$plus = dim(C$plus)[2]
    }
  }


  # Dimension of the tube defining matrix
  rownames(C$minus) <- round(x,3)
  rownames(C$plus) <- round(x,3)

  # Cnames <- list(x = round(x,3), c = 1:dC$minus[2])

  #-----------------------------------------------------------------------------
  # Main code of the algorithm
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
  # Estimate the preimage
  if(Preimage.method$name == "thickening"){
    hatmu1C = PreimageC_thickening(
                    C = C, hatmu = hatmu, hatsigma = hatsigma,
                    tN = tN, kN  = Preimage.method$kN)
  }else if(Preimage.method$name == "true"){
    hatmu1C = PreimageC_thickening(
                    C = C, hatmu = Preimage.method$mu, hatsigma = hatsigma,
                    tN = tN, kN  = Preimage.method$kN)
  }else if(Preimage.method$name == "SCB"){
    hatmu1C = list(minus = rep(TRUE, length(hatmu)),
                   plus  = rep(TRUE, length(hatmu)))
  }else if(Preimage.method$name == "Storey.iid"){
    hatmu1C = Preimage.method$m0
  }

  # Get the modified weights and knots
  new = adapt_knotes_weights(x = x, crit.set = hatmu1C,
                             knots = q.method$knots,
                             I_weights = q.method$I_weights)
  q.method$I_weights <- new$I_weights
  q.method$knots     <- new$knots
  #-----------------------------------------------------------------------------
  # Estimate the quantile for the SCoPES
  fair.q = NULL
  if(q.method$name == "mboot"){
    q = MultiplierBootstrapSplit(alpha   = alpha,
                                 R       = q.method$R,
                                 minus   = hatmu1C$minus,
                                 plus    = hatmu1C$plus,
                                 Mboots  = q.method$Mboots,
                                 method  = q.method$Boottype,
                                 weights = q.method$weights)$q
    if(is.infinite(q)){
      q = 0
    }
  }else if(q.method$name == "fair.mboot"){
    samples = MultiplierBootstrapSplit(alpha   = alpha,
                                       R       = q.method$R,
                                       minus   = hatmu1C$minus,
                                       plus    = hatmu1C$plus,
                                       Mboots  = q.method$Mboots,
                                       method  = q.method$Boottype,
                                       weights = q.method$weights)$samples
    if( !(all(!hatmu1C$minus) && all(!hatmu1C$plus)) ){
      fair.q = fair_quantile_boot(alpha, samples = samples, x = x,
                                  crit.set  = hatmu1C,
                                  knots     = q.method$knots,
                                  I_weights = q.method$I_weights,
                                  type      = q.method$type,
                                  alpha_up  = alpha*(length(q.method$knots) - 1),
                                  maxIter   = q.method$maxIter,
                                  tol       = alpha / 100)
        q = fair.q$u(x)
      }else{
          fair.q = NULL
          q = rep(0, length(x))
      }
  }else if(q.method$name == "gauss.iid"){
    q = maxGauss_quantile(p = 1 - alpha, muC = hatmu1C)
  }else if(q.method$name == "t.iid"){
    q = maxT_quantile(p = 1 - alpha, muC = hatmu1C, df = q.method$df)
  }

  #-----------------------------------------------------------------------------
  # Get the estimates of the excursion sets
  Iminus = dim(C$minus)[2]
  Iplus  = dim(C$minus)[2]

  # Initialize matrizes for the CoPE sets and for the rejections
  hatLC        <- t(t(matrix(FALSE, length(hatmu), Iminus)))
  hatUC        <- t(t(matrix(FALSE, length(hatmu), Iplus)))
  Lcontain_loc <- t(t(matrix(FALSE, length(hatmu), Iminus)))
  Ucontain_loc <- t(t(matrix(FALSE, length(hatmu), Iplus)))
  DetectL  <- t(t(matrix(NA, length(hatmu), Iminus)))
  DetectU  <- t(t(matrix(NA, length(hatmu), Iplus)))

  for(k in 1:max(num.levels$minus, num.levels$plus)){
    # Lower excursion set
    if(k <= num.levels$minus){
      # Compute the lower and upper excursion confidence set, if required
      hatL <- hatLc(C$minus[, k], hatmu, hatsigma, tN, q,
                    mu = mu, inclusion = inclusion$L )
      # Fill the hatLC etc matrix and save where correct inclusions appear
      if(!is.null(mu)){
        hatLC[, k]        <- hatL$hatLc
        Lcontain_loc[, k] <- hatL$contain
        DetectL[, k]      <- hatL$detect
      }else{
        hatLC[, k] <- hatL
      }
    }
    # Upper excursion set
    if(k <= num.levels$plus){
      # Compute the lower and upper excursion confidence set, if required
      hatU <- hatUc(C$plus[, k], hatmu, hatsigma, tN, q,
                    mu = mu, inclusion = inclusion$U )
      # Fill the hatUC matrix and save where correct inclusions appear
      if(!is.null(mu)){
        hatUC[, k]        <- hatU$hatUc
        Ucontain_loc[, k] <- hatU$contain
        DetectU[, k]      <- hatU$detect
      }else{
        hatUC[, k] <- hatU
      }
    }
  }

  if(!is.null(mu)){
    # Set the col and row names for the output variables
    rownames(hatLC) <- rownames(hatUC) <- rownames(C$minus)

    Lcontain <- all(Lcontain_loc)
    Ucontain <- all(Ucontain_loc)
    contain  <- Lcontain && Ucontain

    t_detect <- sum(t(DetectL), na.rm = T) + sum(t(DetectU), na.rm = T)
    f_detect <- sum(t(DetectL) == 0, na.rm = T) + sum(t(DetectU) == 0, na.rm = T)

    evaluation <- list()
    evaluation$contain  = contain
    evaluation$Lcontain = Lcontain
    evaluation$Ucontain = Ucontain
    evaluation$Lcontain_loc = Lcontain_loc
    evaluation$Ucontain_loc = Ucontain_loc
    evaluation$Ldetect  = DetectL
    evaluation$Udetect  = DetectU
    evaluation$num_true_detect  = t_detect
    evaluation$num_false_detect = f_detect

    return(list(hatLC = hatLC, hatUC = hatUC, q = q, hatmu1C = hatmu1C,
                q.method = q.method, Preimage.method = Preimage.method,
                evaluation = evaluation, fair.q = fair.q,
                mu = mu, x = x, tN = tN, hatmu = hatmu, hatsigma = hatsigma,
                C = C))
  }else{
    return(list(hatLC = hatLC, hatUC = hatUC, q = q, hatmu1C = hatmu1C,
                q.method = q.method, Preimage.method = Preimage.method,
                x = x, tN = tN, hatmu = hatmu, hatsigma = hatsigma,
                C = C))
  }
}


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams hatLc
#' @param scopes (list) output of the function SCoPES(...).
#' @param index_C (numeric vector) of length two which contains the index of the
#' lower and the upper excursion set which should be ploted.
#' @param statistic (numeric vector)
#' @param ... Options for plot().
#'  \itemize{
#'   \item L string either "inequal" or "equal" determining whether the lower
#'   excursion set uses < c or <=c.
#'   \item U string either "inequal" or "equal" determining whether the upper
#'   excursion set uses > c or >=c.
#' }
#' @return Plot of the elements within a specified upper and lower
#' excursion set.
#' @export
plot_SCoPES <- function(scopes, index_C = c(1,1),
                        mu = NULL, statistic = NULL, ...){
  # If no statistic is specified use the estimator
  if(is.null(statistic)){
    y = scopes$hatmu
  }else{
    y = statistic
  }
  # Get the correct
  # Get a color vector indicating the "red" (upper excursions) and
  # the blue (lower excursions) set
  colVec <- rep("black", length(scopes$hatmu))
  colVec[scopes$hatLC[, index_C[1]]] <- "blue"
  colVec[scopes$hatUC[, index_C[2]]] <- "red"
  plot(scopes$x, y, col = colVec, ...)
  if( any(scopes$C$minus[, index_C[1]] != scopes$C$plus[, index_C[2]] ) ){
    lines(scopes$x, scopes$C$minus[, index_C[1]], lty = 1, col = "blue")
    lines(scopes$x, scopes$C$plus[, index_C[2]], lty = 1, col = "red")
  }else{
    lines(scopes$x, scopes$C$minus[, index_C[1]], lty = 1, col = "orchid3")
  }
}


#' This functions computes the generalization of the Kac-Rice formula
#' from ...
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

    q = fair_quantile_EEC_t(alpha = alpha,
                            tau   = q.method$tau,
                            x     = range(q.method$knots),
                            df        = q.method$df,
                            crit.set  = crit.set,
                            knots     = q.method$knots,
                            I_weights = q.method$I_weights,
                            type      = "linear",
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
