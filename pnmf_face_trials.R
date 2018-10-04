## ---- Setup ----
library(tidyverse)
library(NMF)
library(pNMF)
#library(hNMF)
## Register the new routines
setNMFMethod("PNMF", pNMF::PNMF)
setNMFMethod("PNMFKL", pNMF::PNMFKL)
setNMFMethod("PNMFO", pNMF::PNMFO)

# projective gradient for comparison
#setNMFMethod("PGNMF", hNMF::PGNMF)

data(faces)
trialrank <- 16
## ---- RandomFaces ----
## Choose a random 16 faces
j <- sample_n(as.data.frame(faces), 16)

ggdf <- NMFimage2df(j, 32, 32)

ggplot(ggdf, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()

## ---- Choices ----
nmfAlgorithm()
nmfSeed()
## ---- Brunet ----
faces.brunet.ica <- nmf(t(faces), rank = trialrank, method='brunet', seed="ica" )
faces.brunet.nndsvd <- nmf(t(faces), rank = trialrank, method='brunet', seed="nndsvd" )

## ---- BrunetRandomSeed ----
faces.brunet.random <- nmf(t(faces), rank = trialrank, method='brunet', seed="random", nrun=20 )

## ---- Offset ----
faces.offset.ica <- nmf(t(faces), rank = trialrank, method='offset', seed="ica" )
faces.offset.nndsvd <- nmf(t(faces), rank = trialrank, method='offset', seed="nndsvd" )

faces.offset.random <- nmf(t(faces), rank = trialrank, method='offset', seed="random", nrun=20 )

## ---- BrunetPlotICA ----
plotNMFBasis <- function(ff) {
  faces.basis.brunet <- NMF::basis(ff)

  ggfbb <- NMFimage2df(t(faces.basis.brunet), 32,32)
  ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
    geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
    coord_fixed()
}
plotNMFBasis(faces.brunet.ica)
## ---- BrunetPlotNNDSVD ----
plotNMFBasis(faces.brunet.nndsvd)

## ---- BrunetPlotRandom ----
plotNMFBasis(faces.brunet.random)

## ---- OffsetPlotICA ----
faces.basis.offset <- NMF::basis(faces.offset.ica)

ggfbb <- NMFimage2df(t(faces.basis.offset), 32,32)
ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()

## ---- OffsetPlotNNDSVD ----
plotNMFBasis(faces.offset.nndsvd)

## ---- OffsetPlotRandom ----
plotNMFBasis(faces.offset.random)

## ---- Projective ----
faces.pnmf.ica <- nmf(t(faces), rank = trialrank, method='PNMF', seed="ica" )
faces.pnmf.nndsvd <- nmf(t(faces), rank = trialrank, method='PNMF', seed="nndsvd" )

faces.pnmf.random <- nmf(t(faces), rank = trialrank, method='PNMF', seed="random", nrun=20 )

## ---- ProjectiveOrtho ----
faces.pnmfo.ica <- nmf(t(faces), rank = trialrank, method='PNMFO', seed="ica" )
faces.pnmfo.nndsvd <- nmf(t(faces), rank = trialrank, method='PNMFO', seed="nndsvd" )

faces.pnmfo.random <- nmf(t(faces), rank = trialrank, method='PNMFO', seed="random", nrun=20 )

## ---- ProjectiveKL ----
faces.pnmfkl.ica <- nmf(t(faces), rank = trialrank, method='PNMFKL', seed="ica" )
faces.pnmfkl.nndsvd <- nmf(t(faces), rank = trialrank, method='PNMFKL', seed="nndsvd" )

faces.pnmfkl.random <- nmf(t(faces), rank = trialrank, method='PNMFKL', seed="random", nrun=20 )


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
                           seed=list(method="nndsvd", densify="average"))
faces.pnmf.nndsvdA <- nmf(t(faces), rank = trialrank,
                         method='PNMF',
                         seed=list(method="nndsvd", densify="average") )

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

## ---- Last ----
