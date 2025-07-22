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
sim_SCoPES <- function(Msim, N, alpha, C, q.method, model, I = NULL,
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
    q.method$R = R

    # Estimate the mean and sd
    hatmu = rowMeans(Y)
    if(model$truesigma){
      hatsigma = model$sigma(model$x)
    }else{
      hatsigma = apply(Y, 1, sd)
    }

    if(!is.null(q.method$mu1Cest)){
      if(q.method$mu1Cest == "m0"){
        pvals <- apply(Y, 1, function(v) t.test(x = v,
                                                alternative = "two.sided",
                                                conf.level = 1-alpha)$p.value)
        q.method$m0 = min(2*sum(pvals >= 0.5), length(pvals))
      }
    }

    # Get the SCoPES for the method
    res_m <- SCoPES(alpha = alpha, C = C, x = model$x, hatmu = hatmu,
                    hatsigma = hatsigma, tN = 1 / sqrt(N),
                    q.method = q.method,
                    inclusion = inclusion, mu = model$mu(model$x))

    # Save the useful variables from the simulation
    hatq[m]    <- res_m$q
    if(!is.null(q.method$mu1Cest)){
      if(q.method$mu1Cest == "m0"){
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



#' This functions simulates SCBs corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_SCBs <- function(Msim, Nvec = c(20, 50, 100, 200),
                     x, alpha = 0.1, q.method,
                     model, mu.model, sd.model = NULL,
                     est_tau = TRUE){

    local.cov    <- global.cov <- list()
    quantile.est <- tau.est    <- list()
    Timing       <- rep(NA, length(Nvec))
    width.L1 <- width.L2 <- NULL

    subI <- sub.intervals(x, q.method$knots,
                          list(minus = rep(TRUE, length(x)),
                               plus  = rep(TRUE, length(x))))$subI

    # Initialize the q.method list for the data.
    q.method.Y <- q.method

    #-------------------------------------------------------------------------------
    # Simulate the
    for(n in 1:length(Nvec)){
      N  = Nvec[n]
      tN = 1 / sqrt(N)

      if(is.null(sd.model)){
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

        if(is.null(sd.model)){
          # Use sample variance as estimate
          sdY = apply(Y, 1, sd)
        }else{
          sdY = sd.model(x = x)
        }

        # Generate residuals
        R = (Y - rowMeans(Y)) / sdY
        q.method.Y$R = R

        # Change the tau function to the estimate from the sample or keep the
        # truth
        if(est_tau){
          q.method.Y$tau = tau_est(R, x)
          tau = q.method.Y$tau
        }else{
          tau = q.method.Y$tau
        }

        # Get the confidence band
        flag = T
        tryCatch(SCB <- fairSCB(alpha, hatmu = mY, hatrho = sdY, tN = 1/sqrt(N),
                               x = x, q.method = q.method.Y, mu = mu.model(x)),
                 error = function(e){flag <<- F})

        #-------------------------------------------------------------------------------
        # Get the coverage of the band
        if(flag){
          global.cov[[n]][m] <- SCB$glob.cov

          local.cov[[n]][,m] <- unlist(lapply(subI, function(l){
            all(SCB$loc.cov[l])
          } ))

        # Get the Euler characteristic of the excursion sets above

          # Save other interesting quantities
          quantile.est[[n]][,m] <- SCB$SCB$q
          tau.est[[n]][,m]      <- tau(x)

          width.L1 <- rbind(width.L1, SCB$width[, 1])
          width.L2 <- rbind(width.L2, SCB$width[, 2])
        }
      }
      Ie <- Sys.time()
      Timing[n] <- Ie - Ia
    }

    #---------------------------------------------------------------------------
    # Create a list with the results
    NI = length(q.method$knots)-1
    na.sims <- vapply(1:length(local.cov), function(l)
                      sum(is.na(local.cov[[l]][1,])),
                      FUN.VALUE = 0.1)
    cov.res <- vapply(1:length(local.cov), function(l)
                      rowMeans(local.cov[[l]], na.rm = TRUE),
                      FUN.VALUE = seq(0, 1, length.out = NI))

    cov.res <- rbind(cov.res, vapply(1:length(local.cov), function(l)
                     mean(global.cov[[l]], na.rm = TRUE), FUN.VALUE = 1))

    cov.res.sd <- vapply(1:length(local.cov), function(l)
                          sqrt(matrixStats::rowVars(1*local.cov[[l]], na.rm = TRUE)),
                          FUN.VALUE = seq(0, 1, length.out = NI))

    cov.res.sd <- rbind(cov.res.sd, vapply(1:length(local.cov), function(l)
                         sd(1*global.cov[[l]], na.rm = TRUE), FUN.VALUE = 1))

    rownames(cov.res) <- rownames(cov.res.sd) <- c(1:NI, "global")
    colnames(cov.res) <- colnames(cov.res.sd) <- Nvec

    return(list(coverage = cov.res,
             coverage.sd = cov.res.sd,
             local.cov   = local.cov,
             na.sims     = na.sims,
             quantiles   = quantile.est,
             tau         = tau.est,
             L1          = width.L1,
             L2          = width.L2,
             time        = Timing,
             simSD       = sqrt(alpha*(1-alpha)/Msim)))
}



#' This functions simulates SCBs corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_SCBs_var <- function(Msim, Nvec = c(20, 50, 100, 200),
                         x, alpha = 0.1, q.method, model, sd.model,
                         est_tau = TRUE){

  local.cov    <- global.cov <- list()
  lb.est       <- ub.est <- tau.est    <- list()
  Timing       <- rep(NA, length(Nvec))

  subI <- sub.intervals(x, q.method$knots,
                        list(minus = rep(TRUE, length(x)),
                             plus  = rep(TRUE, length(x))))$subI

  # Initialize the q.method list for the data.
  q.method.Y <- q.method

  #-------------------------------------------------------------------------------
  # Simulate the
  for(n in 1:length(Nvec)){
    N  = Nvec[n]
    q.method.Y$df = N-1

    local.cov[[n]]  <- matrix(NA, length(subI), Msim)
    global.cov[[n]] <- rep(NA, Msim)

    lb.est[[n]] <- ub.est[[n]] <- matrix(NA, length(x), Msim)
    tau.est[[n]] <- matrix(NA, length(x), Msim)

    Ia <- Sys.time()
    for(m in 1:Msim){
      # Generate a sample
      Y = model(N,  x = x)$values

      # Estimate sd
      sdY = apply(Y, 1, sd)

      # Generate residuals
      R = (Y - rowMeans(Y)) / sdY
      q.method.Y$R = R

      # Change the tau function to the estimate from the sample or keep the
      # truth
      if(est_tau){
        tau = q.method.Y$tau_est(R, x)
      }else{
        tau = q.method.Y$tau
      }
      q.method.Y$tau = tau

      # Get the confidence band
      flag = T
      tryCatch(SCB <- fairSCB_var(alpha, hatvar = sdY^2, x = x,
                                  q.method = q.method.Y, type = "two-sided",
                                  true_var = sd.model(x)^2),
               error = function(e){flag <<- F})

      #-------------------------------------------------------------------------------
      # Get the coverage of the band
      if(flag){
        global.cov[[n]][m] <- SCB$glob.cov

        local.cov[[n]][,m] <- unlist(lapply(subI, function(l){
          all(SCB$loc.cov[l])
        } ))
        # Save other interesting quantities
        lb.est[[n]][,m]  <- SCB$lb(x)
        ub.est[[n]][,m]  <- SCB$ub(x)
        tau.est[[n]][,m] <- tau(x)
      }
    }
    Ie <- Sys.time()
    Timing[n] <- Ie - Ia
  }

  #---------------------------------------------------------------------------
  # Create a list with the results
  NI = length(q.method$knots)-1
  na.sims <- vapply(1:length(local.cov), function(l)
    sum(is.na(local.cov[[l]][1,])),
    FUN.VALUE = 0.1)
  cov.res <- vapply(1:length(local.cov), function(l)
    rowMeans(local.cov[[l]], na.rm = TRUE),
    FUN.VALUE = seq(0, 1, length.out = NI))

  cov.res <- rbind(cov.res, vapply(1:length(local.cov), function(l)
    mean(global.cov[[l]], na.rm = TRUE), FUN.VALUE = 1))

  cov.res.sd <- vapply(1:length(local.cov), function(l)
    sqrt(matrixStats::rowVars(1*local.cov[[l]], na.rm = TRUE)),
    FUN.VALUE = seq(0, 1, length.out = NI))

  cov.res.sd <- rbind(cov.res.sd, vapply(1:length(local.cov), function(l)
    sd(1*global.cov[[l]], na.rm = TRUE), FUN.VALUE = 1))

  rownames(cov.res) <- rownames(cov.res.sd) <- c(1:NI, "global")
  colnames(cov.res) <- colnames(cov.res.sd) <- Nvec

  return(list(coverage = cov.res,
              coverage.sd = cov.res.sd,
              local.cov   = local.cov,
              na.sims     = na.sims,
              lb          = lb.est,
              ub          = ub.est,
              tau         = tau.est,
              time        = Timing,
              simSD       = sqrt(alpha*(1-alpha)/Msim)))
}



#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_LinModelInference <- function(Msim,
                                  NVec = c(20, 50, 100, 200),
                                  x, alpha = 0.1, q.method,
                                  model, sigma.est = TRUE){

  local.cov    <- global.cov <- list()
  quantile.est <- tau.est    <- list()
  Timing       <- rep(NA, length(Nvec))

  subI <- sub.intervals(x, q.method$knots,
                        list(minus = rep(TRUE, length(x)),
                             plus  = rep(TRUE, length(x))))$subI

  # Initialize the q.method list for the data.
  q.method.Y <- q.method

  #-------------------------------------------------------------------------------
  # Simulate the
  for(n in 1:length(Nvec)){
    N  = Nvec[n]
    tN = 1 / sqrt(N)

    if(sigma.est){
      q.method.Y$df = N-1
    }

    local.cov[[n]]  <- matrix(NA, length(subI), Msim)
    global.cov[[n]] <- rep(NA, Msim)

    quantile.est[[n]] <- matrix(NA, length(x), Msim)
    tau.est[[n]]      <- matrix(NA, length(x), Msim)

    # Choose whether sigma is estimated or not
    if(sigma.est){
      sigma = "estimate"
    }else{
      sigma = model$sigma
    }

    Ia <- Sys.time()
    for(m in 1:Msim){
      # Generate the two samples
      Y = model(N,  x = x)$values

      Y1 = sigma * model[[1]]$noise.model(N1, x)$values + model$mu1(x, a1)
      Y2 = sigma * model[[2]]$noise.model(N2, x)$values + model$mu1(x, a1) + model$mu2(x, a2)

      # Estimate the parameters
      linMod_Y = ptw_linModel(Y = cbind(Y1, Y2), X = model$X,
                              sigma = sigma)

      hattau = tau_est(R = linMod_Y$residuals/linMod_Y$sigma, x, df = 40)
      tau    = model$tau

      # Get the confidence band
      flag = T
      tryCatch(SCB <- fairSCB(alpha, hatmu = mY, hatrho = sdY, tN = 1/sqrt(N),
                              x = x, q.method = q.method.Y, mu = mu.model(x)),
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


#' This functions computes the SCoPES corresponding to an estimator and a set
#' of functions given as a matrix with columns being the cut-off functions.
#'
#' @inheritParams SCoPES
#' @return Standard error under the assumption the data is Gaussian
#' @export
sim_SCB_sd <- function(Msim,
                       Nvec = c(20, 50, 100, 200),
                       x, alpha = 0.1, q.method,
                       model, type = "two-sided", estimate.tau = TRUE){

  local.cov    <- global.cov <- list()
  lb.est <- ub.est <- tau.est    <- list()
  Timing       <- rep(NA, length(Nvec))

  subI <- sub.intervals(x, q.method$knots,
                        list(minus = rep(TRUE, length(x)),
                             plus  = rep(TRUE, length(x))))$subI

  # Initialize the q.method list for the data.
  q.method.Y <- q.method

  #-------------------------------------------------------------------------------
  # Simulate the
  for(n in 1:length(Nvec)){
    N  = Nvec[n]

    q.method.Y$df = N-1

    local.cov[[n]]  <- matrix(NA, length(subI), Msim)
    global.cov[[n]] <- rep(NA, Msim)

    tau.est[[n]] <- lb.est[[n]]  <- ub.est[[n]] <- matrix(NA, length(x), Msim)

    Ia <- Sys.time()
    for(m in 1:Msim){
      # Generate a sample
      Y = model$model(N,  x = x)$values

      # Estimate the mean and sd
      mY   = rowMeans(Y)
      varY = matrixStats::rowVars(Y)
      # plot(x, varY)
      # lines(x, model$sd.model(x)^2)

      if(estimate.tau){
        # Generate residuals
        R = (Y - rowMeans(Y)) / sqrt(varY)
        q.method.Y$R = R

        # Estimate tau
        tau = q.method.Y$tau_est(R, x)
      }else{
        tau = model$tau
      }

      # Change the tau function to the estimate from the sample
      q.method.Y$tau = tau

      # Get the confidence band
      flag = T
      tryCatch(SCB <- fairSCB_var(alpha = alpha, hatvar = varY, x = x,
                              q.method = q.method.Y, type = type,
                              true_var = model$sd.model(x)^2),
               error = function(e){flag <<- F})

      # plot(NULL, xlim = range(x), ylim = range(c(SCB$SCB[SCB$SCB!=Inf], model$sd.model(x)^2)),
      #      main= paste("Covering Truth = ", SCB$glob.cov) )
      # lines(x, model$sd.model(x)^2, col = "black")
      # lines(x, SCB$SCB[,"low"], col = "blue", lty = 2)
      # lines(x, SCB$SCB[,"up"], col = "red", lty = 2)

      #-------------------------------------------------------------------------------
      # Get the coverage of the band
      if(flag){
        global.cov[[n]][m] <- SCB$glob.cov

        local.cov[[n]][,m] <- unlist(lapply(subI, function(l){
          all(SCB$loc.cov[l])
        } ))
        # Save other interesting quantities
        lb.est[[n]][,m]  <- SCB$SCB[,"low"]
        ub.est[[n]][,m]  <- SCB$SCB[,"up"]
        tau.est[[n]][,m] <- tau(x)
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
              lb          = lb.est,
              ub          = ub.est,
              tau         = tau.est,
              time        = Timing,
              simSD       = sqrt(alpha*(1-alpha)/Msim)))
}
