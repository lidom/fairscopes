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


