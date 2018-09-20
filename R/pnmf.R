
#' Projective nonnegative matrix factorization
#'
#' @param X Input data matrix
#' @param nmfMod NMF model from the NMF package
#' @param tol tolerance for stopping criteria
#' @param maxIter Maximum number of iterations
#' @param checkDivergence
#'
#' @return Fitted NMF model, as defined in NMF package.
#' @export
#'
#' @importFrom NMF basis coef
#' @examples
#' NMF::register(PNMF)
#'
PNMF <- function (X, nmfMod, tol = 1e-5, maxIter = 500, checkDivergence = TRUE) {
  # Initialization
  startTime <- proc.time()[3]
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
    diffW <- norm(W_old-W, 'fro') / norm(W_old, 'fro')
    if (diffW < tol) {
      break
    }
  }
  gc()
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- t(W)
  return(nmfMod)
}
