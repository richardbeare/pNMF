## ---- Setup ----
library(tidyverse)
library(NMF)
library(pNMF)
#library(hNMF)
## Register the new routines
setNMFMethod("PNMF", pNMF::PNMF)
setNMFMethod("PNMF2", pNMF::PNMF2)
setNMFMethod("PNMFKL", pNMF::PNMFKL)
setNMFMethod("PNMFO", pNMF::PNMFO)
setNMFMethod("PNMFO2", pNMF::PNMFO2)

# projective gradient for comparison
#setNMFMethod("PGNMF", hNMF::PGNMF)

data(faces)
trialrank <- 16

nn <- apply(faces, MARGIN=1, norm, type="2")
faces <- sweep(faces, MARGIN=1, nn, "/")

## parallel options
oo <- list(parallel=10, verbose=TRUE)


## ---- RandomFaces ----
## Choose a random 16 faces
set.seed(233)
j <- sample_n(as.data.frame(faces), 16)

ggdf <- NMFimage2df(j, 32, 32)
rrng <- quantile(ggdf$Brightness, c(0.05, 0.95))

ggplot(ggdf, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) +
  scale_fill_gradient(low="black", high="white", limits=rrng, oob=scales::squish) +
  coord_fixed()

## ---- Choices ----
nmfAlgorithm()
nmfSeed()
## ---- Brunet ----
faces.brunet.ica <- nmf(t(faces), rank = trialrank, method='brunet', seed="ica", nrun=1 )
faces.brunet.nndsvd <- nmf(t(faces), rank = trialrank, method='brunet',
                           seed="nndsvd", nrun=1 )

## ---- BrunetRandomSeed ----
faces.brunet.random <- nmf(t(faces), rank = trialrank, method='brunet', seed="random",
                           nrun=20, .pbackend="mc", .options=oo)

## ---- Offset ----
faces.offset.ica <- nmf(t(faces), rank = trialrank, method='offset', seed="ica", nrun=1 )
faces.offset.nndsvd <- nmf(t(faces), rank = trialrank, method='offset',
                           seed="nndsvd",
                           nrun=1)

faces.offset.random <- nmf(t(faces), rank = trialrank, method='offset', seed="random",
                           nrun=20, .pbackend="mc", .options=oo)

## ---- BrunetPlotICA ----
plotNMFBasis <- function(ff) {
  faces.basis.brunet <- NMF::basis(ff)

  ggfbb <- NMFimage2df(t(faces.basis.brunet), 32,32)
  rrng <- quantile(ggfbb$Brightness, c(0.01, 0.99))

  ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
    geom_tile() + facet_wrap(~ID) +
    scale_fill_gradient(low="black", high="white", limits=rrng, oob=scales::squish) +
    coord_fixed()
}
plotNMFBasis(faces.brunet.ica)
## ---- BrunetPlotNNDSVD ----
plotNMFBasis(faces.brunet.nndsvd)

## ---- BrunetPlotRandom ----
plotNMFBasis(faces.brunet.random)

## ---- OffsetPlotICA ----
plotNMFBasis(faces.brunet.ica)
## ---- OffsetPlotNNDSVD ----
plotNMFBasis(faces.offset.nndsvd)

## ---- OffsetPlotRandom ----
plotNMFBasis(faces.offset.random)

## ---- Projective ----
faces.pnmf.ica <- nmf(t(faces), rank = trialrank, method='PNMF', seed="ica", maxIter=50000, nrun=1 )
faces.pnmf.nndsvd <- nmf(t(faces), rank = trialrank, method='PNMF', seed="nndsvd", maxIter=50000, nrun=1 )

faces.pnmf.random <- nmf(t(faces), rank = trialrank, method='PNMF', seed="random",
                         nrun=20, .pbackend="mc", maxIter=50000, .options=oo )

## ---- ProjectiveOrtho ----
faces.pnmfo.ica <- nmf(t(faces), rank = trialrank, method='PNMFO', seed="ica", maxIter=50000, nrun=1 )
faces.pnmfo.nndsvd <- nmf(t(faces), rank = trialrank, method='PNMFO', seed="nndsvd", maxIter=50000, nrun=1 )

faces.pnmfo.random <- nmf(t(faces), rank = trialrank, method='PNMFO', seed="random",
                          nrun=20, .pbackend="mc", .options=oo, maxIter=50000  )

## ---- ProjectiveKL ----
faces.pnmfkl.ica <- nmf(t(faces), rank = trialrank, method='PNMFKL', seed="ica", maxIter=50000, nrun=1 )
faces.pnmfkl.nndsvd <- nmf(t(faces), rank = trialrank, method='PNMFKL', seed="nndsvd", maxIter=50000, nrun=1 )

faces.pnmfkl.random <- nmf(t(faces), rank = trialrank, method='PNMFKL', seed="random", maxIter=50000,
                           nrun=20, .pbackend="mc", .options=oo)


## ---- ProjectivePlotICA ----
plotNMFBasis(faces.pnmf.ica)

## ---- ProjectivePlotNNDSVD ----
plotNMFBasis(faces.pnmf.nndsvd)

## ---- ProjectivePlotRandom ----
plotNMFBasis(faces.pnmf.random)

## ---- ProjectiveOPlotICA ----
plotNMFBasis(faces.pnmfo.ica)

## ---- ProjectiveOPlotNNDSVD ----
plotNMFBasis(faces.pnmfo.nndsvd)

## ---- ProjectiveOPlotRandom ----
plotNMFBasis(faces.pnmfo.random)

## ---- ProjectiveKLPlotICA ----
plotNMFBasis(faces.pnmfkl.ica)

## ---- ProjectiveKLPlotNNDSVD ----
plotNMFBasis(faces.pnmfkl.nndsvd)

## ---- ProjectiveKLPlotRandom ----
plotNMFBasis(faces.pnmfkl.random)

## ---- PNMFAveSVD ----
faces.pnmfo.nndsvdA <- nmf(t(faces), rank = trialrank,
                           method='PNMFO',
                           seed=list(method="nndsvd", densify="average"),
                           nrun=1)
faces.pnmf.nndsvdA <- nmf(t(faces), rank = trialrank,
                         method='PNMF',
                         seed=list(method="nndsvd", densify="average"),
                         nrun=1)

## ---- PlotPNMFAveSVD ----
plotNMFBasis(faces.pnmfo.nndsvdA)
plotNMFBasis(faces.pnmf.nndsvdA)

## ---- RankSelectPNM ----
ff <- "rsel.Rda"
if (file.exists(ff)) {
  load(ff)
} else {
  oo <- list(parallel=10, verbose=TRUE)
  faces.pnmfo.nndsvdA.rs <- nmf(t(faces), rank = 2:50,
                                method='PNMFO',
                                seed=list(method="nndsvd", densify="average"), nrun=1,
                                .pbackend="mc", .options=oo)
  faces.pnmf.nndsvdA.rs <- nmf(t(faces), rank = 2:50,
                               method='PNMF',
                               seed=list(method="nndsvd", densify="average"),nrun=1,
                               .pbackend="mc", .options=oo)
  save(faces.pnmfo.nndsvdA.rs, faces.pnmf.nndsvdA.rs, file=ff)
}
ff <- "rselrandom.Rda"
if (file.exists(ff)) {
  load(ff)
} else {
  tfacesrand <- randomize(t(faces))

  oo <- list(parallel=10, verbose=TRUE)
  faces.pnmfo.nndsvdA.rs.rand <- nmf(tfacesrand, rank = 2:50,
                                method='PNMFO',
                                seed=list(method="nndsvd", densify="average"), nrun=1,
                                .pbackend="mc", .options=oo)
  faces.pnmf.nndsvdA.rs.rand <- nmf(tfacesrand, rank = 2:50,
                               method='PNMF',
                               seed=list(method="nndsvd", densify="average"),nrun=1,
                               .pbackend="mc", .options=oo)
  save(faces.pnmfo.nndsvdA.rs.rand, faces.pnmf.nndsvdA.rs.rand, file=ff)
}

## ---- PlotRankPNMFO ----
rss.rankplot(faces.pnmfo.nndsvdA.rs, faces.pnmfo.nndsvdA.rs.rand)

## ---- PlotRankPNMF ----
rss.rankplot(faces.pnmf.nndsvdA.rs, faces.pnmf.nndsvdA.rs.rand)

## ---- Reconstruction ----
## using trial rank
compressedFaces <- fitted(faces.pnmf.nndsvdA)
which.faces <- attr(j, "row.names")

ggdf <- NMFimage2df(t(compressedFaces[,which.faces]), 32, 32)

rrng <- quantile(ggdf$Brightness, c(0.05, 0.95))

ggplot(ggdf, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) +
  scale_fill_gradient(low="black", high="white", limits=rrng, oob=scales::squish) +
  coord_fixed()

## ---- PNMFAveSVD ----
## rank of 30 looks OK
faces.pnmfo.nndsvdA <- nmf(t(faces), rank = 30,
                           method='PNMFO',
                           seed=list(method="nndsvd", densify="average"),
                           nrun=1)
faces.pnmf.nndsvdA <- nmf(t(faces), rank = 30,
                          method='PNMF',
                          seed=list(method="nndsvd", densify="average"),
                          nrun=1)

## What about novel data? Don't need it for the structural covariance
compressedFaces <- fitted(faces.pnmf.nndsvdA)
which.faces <- attr(j, "row.names")

ggdf <- NMFimage2df(t(compressedFaces[,which.faces]), 32, 32)

rrng <- quantile(ggdf$Brightness, c(0.05, 0.95))

ggplot(ggdf, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) +
  scale_fill_gradient(low="black", high="white", limits=rrng, oob=scales::squish) +
  coord_fixed()

compressedFaces <- fitted(faces.pnmfo.nndsvdA)
which.faces <- attr(j, "row.names")

ggdf <- NMFimage2df(t(compressedFaces[,which.faces]), 32, 32)

rrng <- quantile(ggdf$Brightness, c(0.05, 0.95))

ggplot(ggdf, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) +
  scale_fill_gradient(low="black", high="white", limits=rrng, oob=scales::squish) +
  coord_fixed()

## ---- NovelData ----

## Projective NMF is special, because H = t(W)*X
## X ~ W t(W) X
## Thus, suppose we look at the coefficients for our random faces
H <- coef(faces.pnmf.nndsvdA)

W <- basis(faces.pnmf.nndsvdA)

aface <- 1964

H[,aface]

t(W) %*% t(faces[1964,,drop=FALSE])

## ---- ParallelRank ----
## The NMF package does lots of stuff about parallel computation to
## make the random seeding workflow faster. This doesn't help
## when the seeding is constant, but we want to check rank.
## The answer is to use mclapply as follows
nmfwrapper <- function(rnk,data, methd, sd) {
  res <- nmf(data, rank=rnk, method=methd, seed=sd, nrun=1)
  return(res)
}
multinmf <- mclapply(2:20, nmfwrapper, data=t(faces), methd='PNMFO',
                     sd=list(method="nndsvd", densify="average"), mc.cores=10)

multinmf.rss <- sapply(multinmf, rss, t(faces))
# Then the same for the randomized data.
## ---- Last ----
