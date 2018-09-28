#' Face data base from the matlab pnmf package
#'
#' Originally from the FERET database, https://www.itl.nist.gov/iad/humanid/feret/.
#'
#' @format Faces are 32x32 images, dataset is 2049 rows,
#' 1024 (=32x32) columns
"faces"



#' ggplot2 dataframe from NMF image data
#'
#' @details transform a typical basis data frame to something that can be
#' displayed as images with ggplot2. Each image is a row and the dimensions
#' need to be passed in.
#' @param nmfimage - data for NMF - each row is an image
#' @param rs - number of rows in each image
#' @param cs - number of cols in each image
#'
#' @return data frame for use with ggplot2
#' @export
#'
#' @examples
#' data(faces)
#' j<-NMFimage2df(faces, 32, 32)
#' ggplot(filter(j, ID <= 16), aes(x=Var2, y=33-Var1, fill=Brightness)) +
#' geom_tile() + facet_wrap(~ID) + scale_fill_gradient(low="black", high="white")
NMFimage2df <- function(nmfimage, rs, cs) {
  ID <- rep(1:nrow(nmfimage))
  nn <- as.data.frame(nmfimage)
  coords <- as.data.frame(expand.grid(1:rs, 1:rs))
  coords$PixelID <- colnames(nn)
  nn$ID <- ID
  nn<- tidyr::gather(nn, key=PixelID, value=Brightness, tidyselect::starts_with("V"))
  nn <- dplyr::left_join(nn, coords, by="PixelID")
  nn
}

#' Plotting for Frigyesi rank selection method
#'
#'
#' @param orig result of a rank selection run on original data
#' @param random result of a rank selection run on randomised data
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#'\dontrun{
#'   mkD <- function(NOISE=TRUE) {
#'   n <- 1000 # features
#'   counts <- c(30, 10, 20, 10, 15, 15) # samples
#'   syntheticNMF(n=n, r=counts, offset = NULL, noise = NOISE,
#'                factors = FALSE, seed = 99)
#'  }
#'  k<-mkD()
#'  V.random <- randomize(k)
#'  aheatmap(V.random)
#'  estim.r2 <- nmf(k, 2:20, nrun=30)
#'  estim.r2.random <- nmf(V.random, 2:20, nrun=30)
#'  rss.rankplot(estim.r2, estim.r2.random)
#'}
rss.rankplot <- function(orig, random) {
  rss.orig <- orig$measures$rss
  rss.random <- random$measures$rss
  df.o <- data.frame(rsschange=-diff(rss.orig), rank=2:length(rss.orig), type="orig")
  df.r <- data.frame(rsschange=-diff(rss.random), rank=2:length(rss.random), type="randomized")
  df<-rbind(df.o, df.r)
  ggplot(df, aes(x=rank, y=rsschange, group=type, colour=type)) + geom_line() + geom_point()
}
