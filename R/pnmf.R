
#' Projective nonnegative matrix factorization
#'
#' @param X Input data matrix
#' @param nmfMod NMF model from the NMF package
#' @param tol tolerance for stopping criteria
#' @param maxIter Maximum number of iterations
#' @param verbose Print status messages
#' @return Fitted NMF model, as defined in NMF package.
#' @export
#'
#' @importFrom NMF basis coef
#' @examples
#' setNMFMethod("PNMF", pNMF::PNMF)
#' mkD <- function(NOISE=TRUE) {
#'   n <- 1000 # rows
#'   counts <- c(30, 10, 20, 10, 15, 15) # samples
#'   syntheticNMF(n=n, r=counts, offset = NULL, noise = NOISE,
#'                factors = FALSE, seed = 99)
#' }
#' k<-mkD()
#' V.random <- randomize(k)
#' estim.r2 <- nmf(k, 2:20, method="PNMF", nrun=30)
#' estim.r2.random <- nmf(V.random, 2:20,  method="PNMF", nrun=30)

PNMF <- function (X, nmfMod, tol = 1e-5, maxIter = 500, verbose=FALSE) {
  # Initialization
  #startTime <- proc.time()[3]
  W <- NMF::basis(nmfMod)
  #H <- NMF::coef(nmfMod)
  err <- rep(0, times = maxIter+1)
  err_diff <- Inf

  # Keep initial W0 for divergence criterion
  W0 <- W
  hasDiverged <- FALSE
  XX <- tcrossprod(X)
  for (iter in 1:maxIter) {
    W_old <- W
    # matlab stuff
    #W = W .* (XX*W) ./ (W*(W'*XX*W) + XX*W*(W'*W));
    #W = W ./ norm(W);
    #diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    W <- W * (XX %*% W)/ (W %*% (crossprod(W, XX)%*%W) + XX %*% W %*% crossprod(W))
    W <- W/norm(W, "2")
    diffW <- norm(W_old-W, 'F') / norm(W_old, 'F')
    if (diffW < tol) {
      if (verbose) {
        cat("Convergence in ", iter, " iterations\n")
      }
      break
    }
  }
  gc()
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- crossprod(W, X)
  return(nmfMod)
}
