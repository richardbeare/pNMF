function() {
  library(NMF)
  library(pNMF)
  data(faces)
  faces <- faces/255

  setNMFMethod("PNMFARD", PNMFARD, overwrite = TRUE)
  setNMFMethod("PNMFARD2", PNMFARD2, overwrite = TRUE)

  NMFTOY <-function(X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
    eps <- .Machine$double.eps
    W <- NMF::basis(nmfMod)
    XX <- tcrossprod(X)
    dsigma <- matrix(0, nrow=ncol(W), ncol=ncol(W))
    checkstep <- 100
    sigma <- colSums(W * W)
    globalW <<- W
    globalXX <<- XX
    globalSigma <<- sigma
    globalX <<- X
    return(nmfMod)
  }
  setNMFMethod("TOY", NMFTOY, overwrite = TRUE)
  k <- nmf(t(faces), rank = 64, method = "TOY", seed = "random", verboseN = TRUE)

  W <- globalW
  sigma <- globalSigma
  XX <- globalXX
  X <- globalX
  dsigma <- diag(sigma)
  library(microbenchmark)
  eps <- .Machine$double.eps

  microbenchmark({W1 <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% (t(W) %*% W)) %*% dsigma + W + eps))}, times=100)
  microbenchmark({W2 <- W * ((XX %*% W) %*% dsigma / ((W %*% (t(W) %*% XX %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))}, times=100)
  microbenchmark({W3 <- W * ((XX %*% W) %*% dsigma / ((W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W)) %*% dsigma + W + eps))}, times=100)

  microbenchmark({XXW <- XX %*% W; W4 <- W * (XXW %*% dsigma / ((W %*% (crossprod(W, XX) %*% W) + XXW %*% crossprod(W)) %*% dsigma + W + eps))}, times=100)

  microbenchmark({XXW <- XX %*% W; W5 <- W * (XXW %*% dsigma / ((W %*% crossprod(W, XXW) + XXW %*% crossprod(W)) %*% dsigma + W + eps))}, times=100)

  all(W1 == W3)
  max(abs(W1 - W3))
  max(abs(W1 - W5))

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

  mat_times_diag <- function(mat, diag) {
    for (i in 1:length(diag)) {
      mat[,i] <- mat[,i] * diag[i]
    }
    return(mat)
  }
  mat_times_diag <- function(mat, diag) {
   mat <- sweep(mat, MARGIN=2,diag, FUN = "*", check.margin = FALSE)
   return(mat)
  }

  W <- WW
  XX <- tcrossprod(t(faces))
  sigma <- colSums(W * W)
library(microbenchmark)
  microbenchmark({W3 <- W * (mat_times_diag(XX %*% W, sigma) / (mat_times_diag(W %*% (crossprod(W, XX) %*% W) + XX %*% W %*% crossprod(W), sigma) + W + eps))}, times=100)
  # No improvement


  ## test equivalence
  k1 <- nmf(t(faces), rank = 64, method = "PNMFARD", seed = "nndsvd", verboseN = TRUE)
  k2 <- nmf(t(faces), rank = 64, method = "PNMFARD2", seed = "nndsvd", verboseN = TRUE)
  k1@extra$sigma - k1@extra$sigma
  k1@extra$wnorm - k1@extra$wnorm

  range(basis(k1) - basis(k2))
  }



#' Automatic rank determination PNMF based on euclidean distance.
#'
#' @details Implementation of "Automatic Rank Determination in Projective Nonnegative Matrix Factorization."
#' Zhirong Yang, Zhanxing Zhu, Erkki Oja. In the 9th International Conference on Latent Variable Analysis
#' and Signal Separation (LVA 2010), pages 514-521, St. Malo, France, 2010
#' Derived from
#' matlab code by Z. Yang, https://sites.google.com/site/zhirongyangcs/ardpnmf
#' @param X Input data matrix
#' @param nmfMod NMF model from the NMF package
#' @param tol tolerance for stopping criteria
#' @param maxIter Maximum number of iterations
#' @param verbose Print status messages
#' @return Fitted NMF model, as defined in NMF package.
#' The "extra" slot contains sigma and wnorm
#' @export
#'
#' @seealso PNMFARD
#'
#' @details sigma and wnorm are stored in the "extra" slot of the nmf object
#' @importFrom NMF basis coef NMFStrategy
#' @examples
#' library(NMF)
#' setNMFMethod("PNMFARD", pNMF::PNMFARD)
#' mkD <- function(NOISE=TRUE) {
#'   n <- 1000 # rows
#'   counts <- c(30, 10, 20, 10, 15, 15) # samples
#'   syntheticNMF(n=n, r=counts, offset = NULL, noise = NOISE,
#'                factors = FALSE, seed = 99)
#' }
#' k<-mkD()
#' estim.r2 <- nmf(k, 16, method="PNMFARD", nrun=1, seed="nndsvd")
#' #wnorm and sigma in the extra slot
#' estim.r2@extra
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
    XXW <- XX %*% W
    W <- W * (XXW %*% dsigma / ((W %*% crossprod(W, XXW) + XXW %*% crossprod(W)) %*% dsigma + W + eps))
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
  wnorm <- apply(W, MARGIN = 2, norm, type = "2")
  stuff <- list(sigma = sigma, wnorm = wnorm)
  nmfMod@extra <- c(nmfMod@extra, stuff)
  return(nmfMod)
}

#' Automatic rank determination PNMF based on euclidean distance - large matrix version.
#'
#' @details Implementation of "Automatic Rank Determination in Projective Nonnegative Matrix Factorization."
#' Zhirong Yang, Zhanxing Zhu, Erkki Oja. In the 9th International Conference on Latent Variable Analysis
#' and Signal Separation (LVA 2010), pages 514-521, St. Malo, France, 2010
#' Derived from
#' matlab code by Z. Yang, https://sites.google.com/site/zhirongyangcs/ardpnmf
#' This one is slower, but uses less RAM.
#' @param X Input data matrix
#' @param nmfMod NMF model from the NMF package
#' @param tol tolerance for stopping criteria
#' @param maxIter Maximum number of iterations
#' @param verbose Print status messages
#' @return Fitted NMF model, as defined in NMF package.
#' The "extra" slot contains sigma and wnorm
#' @export
#'
#' @seealso PNMFARD
#' @details sigma and wnorm are stored in the "extra" slot of the nmf object
#' @importFrom NMF basis coef NMFStrategy
#' @examples
#' library(NMF)
#' setNMFMethod("PNMFARD2", pNMF::PNMFARD2)
#' mkD <- function(NOISE=TRUE) {
#'   n <- 1000 # rows
#'   counts <- c(30, 10, 20, 10, 15, 15) # samples
#'   syntheticNMF(n=n, r=counts, offset = NULL, noise = NOISE,
#'                factors = FALSE, seed = 99)
#' }
#' k<-mkD()
#' estim.r2 <- nmf(k, 16, method="PNMFARD2", nrun=1, seed="nndsvd")
#' #wnorm and sigma in the extra slot
#' estim.r2@extra
PNMFARD2 <-function(X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
  eps <- .Machine$double.eps
  W <- NMF::basis(nmfMod)
  #XX <- tcrossprod(X)
  dsigma <- matrix(0, nrow=ncol(W), ncol=ncol(W))
  checkstep <- 100

  for (iter in 1:maxIter) {
    W_old <- W
    sigma <- colSums(W * W)
    diag(dsigma) <- sigma
    CXW <- crossprod(X,  W)
    XXW <- X %*% CXW

    W <- W * (XXW %*% dsigma /
                ((W %*% (crossprod(W, X) %*% CXW) + XXW %*% crossprod(W)) %*% dsigma + W + eps))
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
  wnorm <- apply(W, MARGIN = 2, norm, type = "2")
  stuff <- list(sigma = sigma, wnorm = wnorm)
  nmfMod@extra <- c(nmfMod@extra, stuff)
  return(nmfMod)
}
