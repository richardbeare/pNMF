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
faces.basis.brunet <- NMF::basis(faces.brunet.ica)

ggfbb <- NMFimage2df(t(faces.basis.brunet), 32,32)
ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()

## ---- BrunetPlotNNDSVD ----
faces.basis.brunet <- NMF::basis(faces.brunet.nndsvd)

ggfbb <- NMFimage2df(t(faces.basis.brunet), 32,32)
ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()

## ---- BrunetPlotRandom ----
faces.basis.brunet <- NMF::basis(faces.brunet.random)

ggfbb <- NMFimage2df(t(faces.basis.brunet), 32,32)
ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()

## ---- OffsetPlotICA ----
faces.basis.offset <- NMF::basis(faces.offset.ica)

ggfbb <- NMFimage2df(t(faces.basis.offset), 32,32)
ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()

## ---- OffsetPlotNNDSVD ----
faces.basis.offset <- NMF::basis(faces.offset.nndsvd)

ggfbb <- NMFimage2df(t(faces.basis.offset), 32,32)
ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()

## ---- OffsetPlotRandom ----
faces.basis.offset <- NMF::basis(faces.offset.random)

ggfbb <- NMFimage2df(t(faces.basis.offset), 32,32)
ggplot(ggfbb, aes(x=Var2, y=33-Var1, fill=Brightness)) +
  geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white") +
  coord_fixed()


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

## ---- Last ----
