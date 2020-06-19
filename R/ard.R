function() {
  library(NMF)
  library(pNMF)
  data(faces)
  faces <- faces/255

  setNMFMethod("PNMFARD", PNMFARD, overwrite = TRUE)


  k <- nmf(t(faces), rank = 64, method = "PNMFARD", seed = "random", verboseN = TRUE)
}


#' @export
PNMFARD <-function(X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
  eps <- .Machine$double.eps
  W <- NMF::basis(nmfMod)
  XX <- tcrossprod(X)
  dsigma <- matrix(0, nrow=ncol(W), ncol=ncol(W))

  checkstep <- 100

  for (iter in 1:maxIter) {
    W_old <- W
    sigma <- colSums(W * W)
    WW <<- W
    Sigma <<- sigma
    diag(dsigma) <- sigma
    break
#browser()
        W <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% (t(W) %*% W)) %*% dsigma + W + eps))
    #W <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))
    #W <- W * ((XX %*% W) %*% dsigma / ((W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))
    W <- W/norm(W, type = "2")
    diffW <- norm(W_old-W, type = "F")/norm(W_old, type = "F")
    if (diffW < tol) {
      if (verboseN) {
        message("Convergence in ", iter, " iterations\n")
      }
      break
    }
    if (verboseN) {
      if (iter %% checkstep == 0) {
        l1 <- paste("Iter ", iter)
        l1 <- paste(l1, paste("Diff=", diffW))
        l1 <- paste(l1, paste("obj=", norm(X - W %*% (t(W) %*% X), type = "F")))
        message(l1)
        message(paste("sigma range", max(sigma), min(sigma)))
      }
    }
  }
  #gc()
  if (verboseN) {
    message(iter, " Iterations used")
  }
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- crossprod(W, X)
#browser()

  return(nmfMod)
}
