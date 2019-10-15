
#' Projective nonnegative matrix factorization based on euclidean distance.
#'
#' @details Implementation of "Linear and Nonlinear Projective Nonnegative Matrix Factorization",
#' Z. Yang and E. Oja, IEEE Transactions on Neural Networks. Derived from
#' matlab code by Z. Yang, https://sites.google.com/site/zhirongyangcs/pnmf.
#' @param X Input data matrix
#' @param nmfMod NMF model from the NMF package
#' @param tol tolerance for stopping criteria
#' @param maxIter Maximum number of iterations
#' @param verbose Print status messages
#' @return Fitted NMF model, as defined in NMF package.
#' @export
#'
#' @importFrom NMF basis coef NMFStrategy
#' @examples
#' library(NMF)
#' setNMFMethod("PNMF", pNMF::PNMF)
#' mkD <- function(NOISE=TRUE) {
#'   n <- 1000 # rows
#'   counts <- c(30, 10, 20, 10, 15, 15) # samples
#'   syntheticNMF(n=n, r=counts, offset = NULL, noise = NOISE,
#'                factors = FALSE, seed = 99)
#' }
#' k<-mkD()
#' estim <- nmf(k, 6, method="PNMF", nrun=1)
#'\dontrun{
#' V.random <- randomize(k)
#' estim.r2 <- nmf(k, 2:20, method="PNMF", nrun=30)
#' estim.r2.random <- nmf(V.random, 2:20,  method="PNMF", nrun=30)
#'}

PNMF <- function (X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
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

  # for testing algebra speed
  #bigW <<- W
  #bigX <<- X
  #bigXX <<- XX
  for (iter in 1:maxIter) {
    W_old <- W
    # matlab stuff
    # Note that the factor of 2 in the paper isn't here
    #W = W .* (XX*W) ./ (W*(W'*XX*W) + XX*W*(W'*W));
    #W = W ./ norm(W);
    #diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    W <- W * (XX %*% W)/ (W %*% crossprod(W, (XX%*%W)) + XX %*% W %*% crossprod(W))
    W <- W/norm(W, "2")
    W[W < zerotol] <- 0
    diffW <- norm(W_old-W, 'F') / norm(W_old, 'F')
    if (diffW < tol) {
      if (verboseN) {
        message("Convergence in ", iter, " iterations\n")
      }
      break
    }
  }
  gc()
  if (verboseN) {
    message(iter, " Iterations used")
  }
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- crossprod(W, X)
  return(nmfMod)
}


## matrix algebra sanity checks
function() {
  load("algebra_tests.Rda")
  W <- bigW
  X <- bigX
  XX <- bigXX
  # matlab code
  # (XX*W) ./ (W*(W'*XX*W) + XX*W*(W'*W))
  # direct translation
  directPort <- (XX %*% W) / (W %*% (t(W) %*% XX %*% W) + XX %*% W %*% (t(W) %*% W))
  myport <- (XX %*% W)/ (W %*% (crossprod(W, XX)%*%W) + XX %*% W %*% crossprod(W))
  max(abs(directPort - myport))
  all(directPort == myport)

  ## timing
  system.time(replicate(100, directPort <- (XX %*% W) / (W %*% (t(W) %*% XX %*% W) + XX %*% W %*% (t(W) %*% W))))
  tW <- t(W)
  system.time(replicate(100, directPort2 <- (XX %*% W) / (W %*% (tW %*% XX %*% W) + XX %*% W %*% (tW %*% W))))

  system.time(replicate(100, myport <- (XX %*% W)/ (W %*% (crossprod(W, XX)%*%W) + XX %*% W %*% crossprod(W))))

  ## Is crossprod faster

  system.time(
    for (i in 1:200) {
      h <- X %*% t(X)
    }
  )
  ## default blas
  ##user  system elapsed
  ##344.136   0.162 344.300
  system.time(replicate(200, h <- X %*% t(X)))
  ## Good - replicate isn't worse
  ## user  system elapsed
  ## 336.552   1.264 337.821

  system.time(replicate(200, h <- tcrossprod(X)))
  ## crossprod a lot faster
  ## user  system elapsed
  ## 227.128   1.556 228.684

  ## same tests with openblas
  ## top shows it cranking along with lots of cores
  # system.time(replicate(200, h <- X %*% t(X)))
  #user  system elapsed
  #75.602  44.937  16.315

  system.time(replicate(200, h <- tcrossprod(X)))
  #user  system elapsed
  #30.965  11.253   7.051
  # crossprod still better

  # XX is computed outside the loop
  system.time(replicate(200, directPort <- (XX %*% W) / (W %*% (t(W) %*% XX %*% W) + XX %*% W %*% (t(W) %*% W))))
  # user  system elapsed
  # 2.378   0.028   2.407

  # very slightly faster
  tW <- t(W)
  system.time(replicate(200, directPort2 <- (XX %*% W) / (W %*% (tW %*% XX %*% W) + XX %*% W %*% (tW %*% W))))
  # user  system elapsed
  # 2.325   0.000   2.326

  ## slightly faster again
  system.time(replicate(200, myport <- (XX %*% W)/ (W %*% (crossprod(W, XX)%*%W) + XX %*% W %*% crossprod(W))))
  # user  system elapsed
  # 2.313   0.000   2.313
  system.time(replicate(200, myport2 <- (XX %*% W)/ (W %*% crossprod(W, (XX%*%W)) + XX %*% W %*% crossprod(W))))
  # user  system elapsed
  # 2.286   0.002   2.288

  max(abs(myport-myport2))

  ## Now for the low memory versions
  system.time(replicate(200, {XtXW <- X %*% crossprod(X, W); myportlm <- XtXW / (tcrossprod(W) %*% XtXW + XtXW %*% crossprod(W))}))
  #user  system elapsed
  #15.827  25.033   5.870
  system.time(replicate(200, {XtXW <- X %*% crossprod(X, W); myportlm <- XtXW / (W %*% crossprod(W, XtXW) + XtXW %*% crossprod(W))}))
  # user  system elapsed
  # 3.591   0.010   3.601
  ## much better

  ## checking ortho versions
  system.time(replicate(200, (XX %*% W)/ (W %*% (crossprod(W, XX)%*%W))))
  #user  system elapsed
  #1.536   0.000   1.537
  system.time(replicate(200, (XX %*% W)/ (W %*% crossprod(W, (XX%*%W))) ))
  # user  system elapsed
  # 1.510   0.000   1.509
  system.time(replicate(200, {
                        XtXW <- X %*% crossprod(X, W)
                        UU <- XtXW/(W %*% crossprod(W, XtXW))
  }))
  #user  system elapsed
  #3.555   0.000   3.556
  ## already has all the tricks in it.
}


#' @describeIn PNMF Projective nonnegative matrix factorization based on euclidean distance.
#' @export
PNMF2 <- function (X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
  # Initialization
  #startTime <- proc.time()[3]
  W <- NMF::basis(nmfMod)
  #H <- NMF::coef(nmfMod)
  err <- rep(0, times = maxIter+1)
  err_diff <- Inf

  # Keep initial W0 for divergence criterion
  W0 <- W
  hasDiverged <- FALSE
  for (iter in 1:maxIter) {
    W_old <- W

    XtXW <- X %*% crossprod(X, W)

    W <- W * XtXW / (W %*% crossprod(W, XtXW) + XtXW %*% crossprod(W))
    W <- W/norm(W, "2")
    W[W < zerotol] <- 0
    diffW <- norm(W_old-W, 'F') / norm(W_old, 'F')
    if (diffW < tol) {
      if (verboseN) {
        message("Convergence in ", iter, " iterations\n")
      }
      break
    }
  }
  gc()
  if (verboseN) {
    message(iter, " Iterations used")
  }

  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- crossprod(W, X)
  return(nmfMod)
}

#' Projective orthonormal nonnegative matrix factorization based on euclidean distance.
#'
#' @details Implementation of "Linear and Nonlinear Projective Nonnegative Matrix Factorization",
#' Z. Yang and E. Oja, IEEE Transactions on Neural Networks. Derived from
#' matlab code by Z. Yang, https://sites.google.com/site/zhirongyangcs/pnmf.
#'
#' The PNMFO2 version uses a different ordering of matrix operations that is slower, as
#' more happens in the loop, but reduces the maximum matrix size. No need for a
#' features x features matrix
#'
#' @aliases PNMFO2
#' @param X Input data matrix
#' @param nmfMod NMF model from the NMF package
#' @param tol tolerance for stopping criteria
#' @param maxIter Maximum number of iterations
#' @param verbose Print status messages
#' @return Fitted NMF model, as defined in NMF package.
#' @export
#'
#' @importFrom NMF basis coef NMFStrategy
#' @examples
#' library(NMF)
#' setNMFMethod("PNMFO", pNMF::PNMFO)
#' mkD <- function(NOISE=TRUE) {
#'   n <- 1000 # rows
#'   counts <- c(30, 10, 20, 10, 15, 15) # samples
#'   syntheticNMF(n=n, r=counts, offset = NULL, noise = NOISE,
#'                factors = FALSE, seed = 99)
#' }
#' k<-mkD()
#' estim <- nmf(k, 6, method="PNMF", nrun=1)
#'\dontrun{
#' V.random <- randomize(k)
#' estim.r2 <- nmf(k, 2:20, method="PNMF", nrun=30)
#' estim.r2.random <- nmf(V.random, 2:20,  method="PNMF", nrun=30)
#'}

PNMFO <- function (X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
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

    W <- W * (XX %*% W)/(W %*% crossprod(W, (XX%*%W)))
    W <- W/norm(W, "2")
    W[W < zerotol] <- 0
    diffW <- norm(W_old-W, 'F') / norm(W_old, 'F')
    if (diffW < tol) {
      if (verboseN) {
        message("Convergence in ", iter, " iterations\n")
      }
      break
    }
  }
  if (verboseN) {
    message(iter, " Iterations used")
  }
  gc()
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- crossprod(W, X)
  return(nmfMod)
}


#' @describeIn PNMFO Projective orthonormal nonnegative matrix factorization based on euclidean distance.
#' @export
PNMFO2 <- function (X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
  # Initialization
  #startTime <- proc.time()[3]
  W <- NMF::basis(nmfMod)
  #H <- NMF::coef(nmfMod)
  err <- rep(0, times = maxIter+1)
  err_diff <- Inf

  # Keep initial W0 for divergence criterion
  W0 <- W
  hasDiverged <- FALSE

  for (iter in 1:maxIter) {
    W_old <- W

    XtXW <- X %*% crossprod(X, W)
    UU <- XtXW/(W %*% crossprod(W, XtXW))
    W <- W * UU
    W <- W/norm(W, "2")
    W[W < zerotol] <- 0
    diffW <- norm(W_old-W, 'F') / norm(W_old, 'F')
    if (diffW < tol) {
      if (verboseN) {
        message("Convergence in ", iter, " iterations\n")
      }
      break
    }
  }
  if (verboseN) {
    message(iter, " Iterations used")
  }
  gc()
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- crossprod(W, X)
  return(nmfMod)
}

#' Projective nonnegative matrix factorization based on I-divergence (non-nomalized KL-divergence)
#'
#' @param X Input data matrix
#' @param nmfMod NMF model from the NMF package
#' @param tol tolerance for stopping criteria
#' @param maxIter Maximum number of iterations
#' @param verbose Print status messages
#' @return Fitted NMF model, as defined in NMF package.
#' @export
#'
#' @importFrom NMF basis coef NMFStrategy

#' @examples
#' library(NMF)
#' NMF::setNMFMethod("PNMFKL", pNMF::PNMFKL)
#' mkD <- function(NOISE=TRUE) {
#'   n <- 1000 # rows
#'   counts <- c(30, 10, 20, 10, 15, 15) # samples
#'   syntheticNMF(n=n, r=counts, offset = NULL, noise = NOISE,
#'                factors = FALSE, seed = 99)
#' }
#' k<-mkD()
#' estim <- nmf(k, 6, method="PNMFKL", nrun=1)
#'
#'\dontrun{
#' V.random <- randomize(k)
#' estim.r2 <- nmf(k, 2:20, method="PNMFKL", nrun=30)
#' estim.r2.random <- nmf(V.random, 2:20,  method="PNMF", nrun=30)
#'}
#'
PNMFKL <- function(X, nmfMod, tol = 1e-5, maxIter = 5000, verboseN=FALSE, zerotol=1e-10) {
  W <- NMF::basis(nmfMod)
  Xsum <- rowSums(X)
  dim(Xsum) <- c(nrow(X), 1)
  for (iter in 1:maxIter) {
    W_old <- W
    # matlab version
    #Z = X ./ (W*(W'*X));
    #W = W .* sqrt((Z*(X'*W) +X*(Z'*W)) ...
    #                            ./ bsxfun(@plus, Xsum'*W, bsxfun(@times, Xsum, sum(W))));

    #diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    Z <- X/((W %*% crossprod(W, X)))
    ## slightly messy in order to implement bsxfun without extra stuff
    ## bsxfun(@times, Xsum, sum(W))
    sW <- colSums(W)
    dim(sW) <- c(1, ncol(W))
    bsxfuntimes <- Xsum %*% sW
    XW <- as.vector(crossprod(Xsum, W))
    bsxfunplus <- sweep(bsxfuntimes, MARGIN=2, STATS=XW, FUN="+")
    W <- W * sqrt((Z %*% crossprod(X, W) + X %*% crossprod(Z, W)))/bsxfunplus
    W[W < zerotol] <- 0
    diffW <- norm(W_old - W, "F")/norm(W_old, "F")
    if (diffW<tol) {
      if (verboseN) {
        message("Converged after ", iter , " steps.\n")
      }
      break;
    }
  }
  if (verboseN) {
    message(iter, " Iterations used")
  }
  gc()
  NMF::basis(nmfMod) <- W
  NMF::coef(nmfMod) <- crossprod(W, X)
  return(nmfMod)
}
