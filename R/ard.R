function() {
  library(NMF)
  library(pNMF)
  data(faces)
  faces <- faces/255

  setNMFMethod("PNMFARD", PNMFARD, overwrite = TRUE)


  k <- nmf(t(faces), rank = 64, method = "PNMFARD", seed = "random", verboseN = TRUE)

  W <- WW
  dsigma <- diag(Sigma)
  library(microbenchmark)
  eps <- .Machine$double.eps
  XX <- tcrossprod(t(faces))
  microbenchmark({W1 <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% (t(W) %*% W)) %*% dsigma + W + eps))}, times=100)
  microbenchmark({W2 <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))}, times=100)
  microbenchmark({W3 <- W * ((XX %*% W) %*% dsigma / ((W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))}, times=100)

  all(W1 == W3)
  max(abs(W1 - W3))
  # slower!
  #library(spam)
  #dsigma2 <- diag.spam(Sigma)
  #microbenchmark({W1 <- W * ((XX %*% W) %*% dsigma2 / ((W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma2 + W + eps))}, times=100)

  # create a special version of dsigma that gives the same answer via element-wise multiplication
  dsigma3 <- matrix(rep(Sigma, rep(nrow(XX), length(Sigma))), ncol=length(Sigma))
  all( ((XX %*% W) %*% dsigma) == ((XX %*% W) * dsigma3))
  microbenchmark({W1 <- W * ((XX %*% W) * dsigma3 / ((W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W)) * dsigma3 + W + eps))}, times=100)
  # No difference, and that doesn't include generating sigma3
  library(Matrix)
  dsigma3 <- Diagonal(x=Sigma)
  microbenchmark({W3 <- W * ((XX %*% W) %*% dsigma3 / ((W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma3 + W + eps))}, times=100)
  XXM <- Matrix(XX)
  WM <- Matrix(W)
  microbenchmark({W3 <- WM * ((XXM %*% WM) %*% dsigma3 / ((WM %*% (crossprod(WM, XXM) %*% WM) + XXM %*% WM %*% crossprod(WM)) %*% dsigma3 + WM + eps))}, times=100)
  # Nope - slower. No point messing with the sigma assuming that it stays relatively small
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
    diag(dsigma) <- sigma
    # W <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% (t(W) %*% W)) %*% dsigma + W + eps))
    #W <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))
    W <- W * ((XX %*% W) %*% dsigma / ((W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))
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
  stuff <- list(sigma = sigma, wnorm = )
  nmfMod@extra <- c(nmfMod@extra, stuff)
  return(nmfMod)
}
