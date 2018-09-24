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
