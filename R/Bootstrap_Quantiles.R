#------------------------------------------------------------------------------#
#                                                                              #
#     Bootstrap functions to obtain the quantiles of the maximum of random     #
#     fields                                                                   #
#                                                                              #
#------------------------------------------------------------------------------#
# Required packages:
#      - matrixStats
#
# Contained functions:
#      - MultiplierBootstrapSplit (tested)
#
#------------------------------------------------------------------------------#
# Developer notes:
# - Style guideline included
#------------------------------------------------------------------------------#
#' Multiplier bootstrap estimator for the quantile of the maximum of a random field. The
#' function contains a multiplier-t version based on Telschow Schwartzmann (2019)
#' "Simultaneous confidence bands for non-parametric regression with functional data"
#' and a regular version based on Cherno. For large sample sizes the versions do agree,
#' however, for small sample sizes the bootstrap-t has better covering rates, but a slightly
#' higher variability in the estimate of the quantile than the simple version based on Cherno.
#'
#' @param R array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain.
#' @param Q array of dimension K_1 x ... x K_D x N containing N-realizations of
#'  a random field over a D-dimensional domain. Default NULL, i.e. one sample
#'  case.
#' @param alpha numeric the targeted upper quantile of the maxiumum of the
#'   absolute value of the random field. Default is 0.95.
#' @param Mboots numeric amount of bootstrap replicates. Default is 5e3.
#' @param method string specifies the bootstrap version. Options are "t" and
#'  "regular". Default is "t".
#' @param weights string specifying the multipliers. Options are "gaussian",
#'  "rademacher" and "mammen". Default is "rademacher".
#' @return list with elements
#'  \itemize{
#'   \item z Vector containing the realizations of the bootstrap of the approximation of the maximum of the random field
#'   \item q alpha-quantile value of z
#' }
#' @export
MultiplierBootstrapSplit <- function(R,
                                     plus,
                                     minus,
                                     Q = NULL,
                                     alpha   = 0.05,
                                     Mboots  = 5e3,
                                     method  = "t",
                                     weights = "rademacher" ){

  #---- Check user input and put default values
  # Check R
  if( !is.array( R ) ){
    stop("'R' must be an array.")
  }
  # Check Q
  if( !( is.array( Q ) | is.null( Q ) ) ){
    stop("'Q' must be an array.")
  }

  #---- Check params and put defaults, if missing
  # Check Mboots
  if( is.numeric( Mboots ) ){
    if( Mboots %% 1 != 0 & Mboots <= 0 ){
      stop( "The input 'Mboots' needs to be a positiv natural number." )
    }
  }else{
    stop("The input 'Mboots' needs to be a positiv natural number.")
  }

  # Check alpha
  if( is.numeric( alpha ) ){
    if( alpha <= 0 | alpha >= 1 ){
      stop("The input 'alpha' needs to be a real number between 0 and 1.")
    }
  }else{
    stop("The input 'alpha' needs to be a real number between 0 and 1.")
  }

  # Check method
  if( is.character( method ) ){
    if( !( method %in% c( "t", "regular" ) ) ){
      stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'method'. Options are 't' and 'regular'.")
  }

  # Check weights
  if( is.character( weights ) ){
    if( !( weights %in% c( "gaussian", "rademacher", "mammen" ) ) ){
      stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
    }
  }else{
    stop( "Please, specify a valid choice for 'weights'. Options are 'gaussian', 'rademacher' and 'mammen'.")
  }

  #----- One sample case
    #----- Precompute useful constants
    # dimension of input
    dimR = dim( R )

    # number of samples
    N    = dimR[ length( dimR ) ]
    D    = length( dimR ) - 1

    #----- Simulate multiplier weights
    if( weights == "gaussian" ){
      multiplier <- matrix( rnorm( N * Mboots ), N, Mboots )
    }else if( weights == "rademacher" ){
      multiplier <- matrix( sample( c( -1, 1 ), N * Mboots, replace = T ), N, Mboots )
    }else{
      multiplier <- matrix( sqrt( 5 ) *
                              rbinom( N * Mboots,
                                      1,
                                      ( sqrt( 5 ) - 1 ) / 2 / sqrt( 5 ) ) +
                              ( 1 - sqrt( 5 ) ) / 2,
                            N,
                            Mboots )
    }

      # Compute bootstrap means
      bootMeans <- R %*% multiplier / N

      # Estimate the variance from the sample
      if( method == "regular" ){
        data.sigma <- t(t(sqrt( matrixStats::rowVars( R ) )))
      }else if( method == "t" ){
        bootSecMoments <- R^2 %*% multiplier^2 / N
        # We put an abs here to make sure that no NaNs are produced due to machine precision error.
        data.sigma <- sqrt( ( N / ( N - 1 ) ) * abs( bootSecMoments - bootMeans^2 ) )
      }

      # Compute bootstrap distribution of the maximum
      if(sum(plus) == 1){
        distVec1 <- sqrt( N ) * bootMeans[plus,] / data.sigma[plus,]
      }else if(sum(plus) > 1){
        distVec1 <- sqrt( N ) * apply( t(t(bootMeans[plus,])) / data.sigma[plus,], 2, max )
      }else{
        distVec1 <- rep(-Inf, Mboots)
      }
      if(sum(minus) == 1){
        distVec2 <- -sqrt( N ) * bootMeans[minus,] / data.sigma[minus,]
      }else if(sum(minus) > 1){
        distVec2 <- sqrt( N ) * apply( -t(t(bootMeans[minus,])) / data.sigma[minus,], 2, max )
      }else{
        distVec2 <- rep(-Inf, Mboots)
      }

      if( method == "regular" ){
        data.sigma = as.vector(data.sigma)
      }

      distVec  <- apply( rbind(distVec1, distVec2), 2, max )

  #----- Return quantile and bootstrap distribution
  return( list( z       = distVec,
                q       = quantile( distVec, 1 - alpha, type = 8 ),
                samples = bootMeans / data.sigma * sqrt( N ) ) )
}
