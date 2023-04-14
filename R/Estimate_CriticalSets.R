#------------------------------------------------------------------------------#
#                                                                              #
#     Functions to estimate critical sets for SCoPE sets                       #
#                                                                              #
#------------------------------------------------------------------------------#
# Contained functions:
# - PreimageC_thickening()
#
#------------------------------------------------------------------------------#
# Developer notes:
#
#------------------------------------------------------------------------------#
#' Estimates the preimage of the set indicated by the matrix C using
#' the thickening estimator proposed in Telschow et al (2022)
#' "SCoPES: A versatile framework of simultaneous inference".
#'
#' @param C (numeric matrix) T x nlvl matrix containing nlvl functions for which the preimage
#'          should be computed. At the moment only the first and last function
#'          matter, i.e., they define a tube.
#' @param hatmu (numeric vector) estimate of the target function mu.
#' @param hatsigma (numeric) scaling used to estimate the set mu1_C. Usually it is
#' the pointwise standard deviation of hatmu.
#' @param tN (numeric) is inverse of the rate of the fCLT, i.e., usually 1/sqrt(N)
#' @param kN (numeric) is the thickening of the set (Default log(tN^-2) / 5)
#' @param method (string) "two-sided" or "one-sided" neighborhood around the
#' thresholding functions defined in C. Default "two-sided".
#' @return list with elements
#'  \itemize{
#'   \item "minus" boolean vector indicating which locations are in the set
#'         mu1_C^-
#'   \item "plus" boolean vector indicating which locations are in the set
#'         mu1_C^+
#' }
#' @export
PreimageC_thickening <- function(C, hatmu, hatsigma, tN,
                                 kN  = log(tN^-2) / 5, method ="two-sided"){
  if(method == "two-sided"){
    Sminus = abs(C$minus - hatmu) <= tN * kN * hatsigma
    Splus  = abs(C$plus - hatmu) <= tN * kN * hatsigma
  }else{
    Sminus = (hatmu - C$minus <= tN * kN * hatsigma) & (0 <= hatmu - C$minus)
    Splus  = (C$plus - hatmu <= tN * kN * hatsigma)  & (0 <=  C$plus - hatmu)
  }

  hatmu1_C = list(minus = Sminus, plus  = Splus)

  # Return the estimate of the preimage of C
  return(hatmu1_C)
}
