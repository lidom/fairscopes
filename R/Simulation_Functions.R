#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to Reproduce Simulation Results from the Article
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
# - sim_SCoPES()
# - sim_SCBs()
#------------------------------------------------------------------------------#
# Developer notes:
# - Fix the documentation and man pages
#
#------------------------------------------------------------------------------#
#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_SCoPES <- function(Msim, N, alpha, C, q_method, model, I = NULL,
                       inclusion = list(L = "inequal", U = "inequal")){
  # Dimension of the tube defining matrix
  dC    = dim(C)
  if(is.null(dC)){
    dC <- c(length(C), 1)
  }
  # Get standard values for the indices of C^\pm if not specified
  if(is.null(I)){
    if(SCoPEStype %in% c( "extraction", "relevance", "lrelevance",
                          "equivalence", "lequivalence")){
      Iminus <- seq(1, dC[2], 2)
      Iplus  <- seq(2, dC[2], 2)
    }else{
      Iminus <- 1:dC[2]
      Iplus  <- 1:dC[2]
    }
    I = list(minus = Iminus, plus = Iplus)
  }else if(is.list(I)){
    if(length(I) == 2){
      Iminus = I$minus
      Iminus = I$plus
    }else{
      stop("I must be a list with two entries indicating which columns of
           C belong to C^- and which to C^+")
    }
  }else{
    stop("I must be a list with two entries indicating which columns of
           C belong to C^- and which to C^+")
  }

  contain <- 0
  hatq    <- rep(NaN, Msim)
  hatmu1C <- list(minus = matrix(FALSE, length(model$x), Msim),
                  plus  = matrix(FALSE, length(model$x), Msim))
  Lcontain <- Ucontain <- matrix(FALSE, length(model$x), Msim)
  detectL  <- array(FALSE, dim = c(dC[1], length(Iminus), Msim))
  detectU  <- array(FALSE, dim = c(dC[1], length(Iplus), Msim))
  NDetect_BH <- NDetect_hommel <- NDetect_sidak <- NDetect  <-  matrix(NaN, 2, Msim)

  for(m in 1:Msim){
    # Get the sample for the simulation run
    Y = SampleFields::SignalPlusNoise(N,
                                      x = model$x,
                                      mu = model$mu,
                                      noise = model$noise,
                                      sigma = model$sigma)
    Y = Y$values
    # Get the residuals
    R = Y - rowMeans(Y)
    q_method$R = R

    # Estimate the mean and sd
    hatmu = rowMeans(Y)
    if(model$truesigma){
      hatsigma = model$sigma(model$x)
    }else{
      hatsigma = apply(Y, 1, sd)
    }

    if(!is.null(q_method$mu1Cest)){
      if(q_method$mu1Cest == "m0"){
        pvals <- apply(Y, 1, function(v) t.test(x = v,
                                                alternative = "two.sided",
                                                conf.level = 1-alpha)$p.value)
        q_method$m0 = min(2*sum(pvals >= 0.5), length(pvals))
      }
    }

    # Get the SCoPES for the method
    res_m <- SCoPES(alpha = alpha, C = C, x = model$x, hatmu = hatmu,
                    hatsigma = hatsigma, tN = 1 / sqrt(N),
                    method = q_method,
                    inclusion = inclusion, mu = model$mu(model$x))

    # Save the useful variables from the simulation
    hatq[m]    <- res_m$q
    if(!is.null(q_method$mu1Cest)){
      if(q_method$mu1Cest == "m0"){
        hatmu1C$minus[, m] <- rep(T, length(model$x))
        hatmu1C$plus[, m]  <- rep(T, length(model$x))
      }else{
        hatmu1C$minus[, m] <- res_m$hatmu1C$minus
        hatmu1C$plus[, m]  <- res_m$hatmu1C$plus
      }
    }else{
      hatmu1C$minus[, m] <- res_m$hatmu1C$minus
      hatmu1C$plus[, m]  <- res_m$hatmu1C$plus
    }
    detectL[,, m]      <- res_m$hatLC
    detectU[,, m]      <- res_m$hatUC
    Lcontain[, m]      <- res_m$Lcontain_loc
    Ucontain[, m]      <- res_m$Ucontain_loc
    NDetect[1, m]      <- res_m$NtrueDetect
    NDetect[2, m]      <- res_m$NfalseDetect
    contain <- contain + res_m$contain / Msim

    if(SCoPEStype == "classical"){
      I0 = model$mu(model$x) == C
      I1 = model$mu(model$x) != C
      #---------------------------------------------------------------------------
      # Testing alternatives
      pvals <- apply(Y, 1, function(v) t.test(x = v,
                                              alternative = "two.sided",
                                              conf.level = 1-alpha)$p.value)

      # pvals_adjust = p.adjust(pvals, method = "sidak")
      detect_hommel  = MPT_Detect(alpha, pvals, "hommel")
      detect_sidak = FWE_control(alpha, pvals, "sidak")
      detect_BH    = FDR_control(alpha, pvals, "BH")

      NDetect_hommel[1, m]  <- sum(detect_hommel[I1])
      NDetect_hommel[2, m]  <- sum(detect_hommel[I0])
      NDetect_sidak[1, m] <- sum(detect_sidak[I1])
      NDetect_sidak[2, m] <- sum(detect_sidak[I0])
      NDetect_BH[1, m] <- sum(detect_BH[I1])
      NDetect_BH[2, m] <- sum(detect_BH[I0])
    }
  }

  if(SCoPEStype == "classical"){
    return(list(coverage  = contain, Lcoverage = Lcontain,
                Ucoverage = Ucontain, q = hatq, mu1C = hatmu1C,
                detectL   = detectL, detectU = detectU,
                NDetect       = NDetect,
                NDetect_hommel  = NDetect_hommel,
                NDetect_sidak = NDetect_sidak,
                NDetect_BH    = NDetect_BH))
  }else{
    return(list(coverage = contain, Lcoverage = Lcontain,
                Ucoverage   = Ucontain, q = hatq, mu1C = hatmu1C,
                detectL     = detectL, detectU = detectU, NDetect = NDetect))
  }
}



#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_SCBs <- function(Msim, NVec = c(20, 50, 100, 200),
                     alpha = 0.1, q_method, model, I = NULL){

    local.cov  <- list()
    global.cov <- list()
    Timing     <- rep(NA, length(Nvec))

    #-------------------------------------------------------------------------------
    # Simulate the
    for(n in 1:length(Nvec)){
      N  = Nvec[n]
      tN = 1 / sqrt(N)

      local.cov[[n]]  <- matrix(NA, Nknots, Msim)
      global.cov[[n]] <- rep(NA, Msim)

      Ia <- Sys.time()
      for(m in 1:Msim){
        # Generate a sample and the residuals
        Y = SignalPlusNoise(N, x = x,
                            mu = mu_model,
                            noise = noise_model,
                            sigma = sigma_model)
        Y = Y$values
        R = (Y - rowMeans(Y)) / sigma_model(x)

        # Estimate the mean and sd
        mY  = rowMeans(Y)
        sdY = apply(Y, 1, sd)

        # Estimate tau
        tau = tau_est(R, x)

        # Get the confidence band
        u = fairEEC_z(alpha, knots, tau = tau)

        #-------------------------------------------------------------------------------
        # Get the coverage of the band
        cov_locs = abs(mu_model(x) - mY) <= u$u(x)*sigma_model(x)*tN
        global.cov[[n]][m] <- all(cov_locs)

        local.cov[[n]][,m] <- unlist(lapply(subI, function(l){
          all(cov_locs[l])
        } ))

      }
      Ie <- Sys.time()
      Timing[n] <- Ie - Ia

    }


    loc.cov.res <- vapply(1:length(local.cov), function(l)
      rowMeans(local.cov[[l]]), FUN.VALUE = seq(0, 1, length.out = Nknots))

    loc.cov.res <- rbind(loc.cov.res, vapply(1:length(local.cov), function(l)
      mean(global.cov[[l]]), FUN.VALUE = 1))

    rownames(loc.cov.res) <- c(1:Nknots, "global")
    colnames(loc.cov.res) <- Nvec
}
