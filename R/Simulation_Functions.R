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
                NDetect        = NDetect,
                NDetect_hommel = NDetect_hommel,
                NDetect_sidak  = NDetect_sidak,
                NDetect_BH     = NDetect_BH))
  }else{
    return(list(coverage  = contain, Lcoverage = Lcontain,
                Ucoverage = Ucontain, q = hatq, mu1C = hatmu1C,
                detectL   = detectL, detectU = detectU, NDetect = NDetect))
  }
}



#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_SCBs <- function(Msim, NVec = c(20, 50, 100, 200),
                     x, alpha = 0.1, q_method, model, mu_model, sd_model = NULL){

    local.cov  <- global.cov <- list()
    quantile.est <- tau.est <- list()
    Timing     <- rep(NA, length(Nvec))

    subI <- sub.intervals(x, q_method$knots,
                          list(minus = rep(TRUE, length(x)),
                               plus  = rep(TRUE, length(x))))$subI

    # Initialize the q.method list for the data.
    q.method.Y <- q.method

    #-------------------------------------------------------------------------------
    # Simulate the
    for(n in 1:length(Nvec)){
      N  = Nvec[n]
      tN = 1 / sqrt(N)

      if(is.null(sd_model)){
        q.method.Y$df = N-1
      }

      local.cov[[n]]  <- matrix(NA, length(subI), Msim)
      global.cov[[n]] <- rep(NA, Msim)

      quantile.est[[n]] <- matrix(NA, length(x), Msim)
      tau.est[[n]]      <- matrix(NA, length(x), Msim)

      Ia <- Sys.time()
      for(m in 1:Msim){
        # Generate a sample
        Y = model(N,  x = x)$values

        # Estimate the mean and sd
        mY  = rowMeans(Y)

        if(is.null(sd_model)){
          # Use sample variance as estimate
          sdY = apply(Y, 1, sd)
        }else{
          sdY = sd_model(x = x)
        }

        # Generate residuals

        R = (Y - rowMeans(Y)) / sdY

        # Estimate tau
        tau = q.method$tau.est(R, x)

        # Change the tau function to the estimate from the sample
        q.method.Y$tau.est = tau

        # Get the confidence band
        flag = T
        tryCatch(SCB <- fairSCB(alpha, hatmu = mY, hatrho = sdY, tN = 1/sqrt(N),
                               x = x, q.method = q.method.Y, mu = mu_model(x)),
                 error = function(e){flag <<- F})

        #-------------------------------------------------------------------------------
        # Get the coverage of the band
        if(flag){
          global.cov[[n]][m] <- SCB$glob.cov

          local.cov[[n]][,m] <- unlist(lapply(subI, function(l){
            all(SCB$loc.cov[l])
          } ))
          # Save other interesting quantities
          quantile.est[[n]][,m] <- SCB$q
          tau.est[[n]][,m]      <- tau(x)
        }
      }
      Ie <- Sys.time()
      Timing[n] <- Ie - Ia
    }

    #---------------------------------------------------------------------------
    # Create a list with the results
    na.sims <- vapply(1:length(local.cov), function(l)
                      sum(is.na(local.cov[[l]][1,])),
                      FUN.VALUE = 0.1)
    cov.res <- vapply(1:length(local.cov), function(l)
                      rowMeans(local.cov[[l]], na.rm = TRUE),
                      FUN.VALUE = seq(0, 1, length.out = q.method$Nknots))

    cov.res <- rbind(cov.res, vapply(1:length(local.cov), function(l)
                     mean(global.cov[[l]], na.rm = TRUE), FUN.VALUE = 1))

    cov.res.sd <- vapply(1:length(local.cov), function(l)
                          sqrt(matrixStats::rowVars(1*local.cov[[l]], na.rm = TRUE)),
                          FUN.VALUE = seq(0, 1, length.out = q.method$Nknots))

    cov.res.sd <- rbind(cov.res.sd, vapply(1:length(local.cov), function(l)
                         sd(1*global.cov[[l]], na.rm = TRUE), FUN.VALUE = 1))

    rownames(cov.res) <- rownames(cov.res.sd) <- c(1:q.method$Nknots, "global")
    colnames(cov.res) <- colnames(cov.res.sd) <- Nvec

    return(list(coverage = cov.res,
             coverage.sd = cov.res.sd,
             local.cov   = local.cov,
             na.sims     = na.sims,
             quantiles   = quantile.est,
             tau         = tau.est,
             time        = Timing,
             simSD       = sqrt(alpha*(1-alpha)/Msim)))
}
