
figure <- function(title, p){
  RESULT_FOLDER <- "results"
  ggsave(file.path(RESULT_FOLDER, paste0(make.names(title),".pdf")))
  p <- p + ggtitle(title)
  ggsave(file.path(RESULT_FOLDER, paste0(make.names(title),".png")), dpi=600)
  p
}


figure1 <- function(title_, p, sub_title=FALSE, title=!sub_title){
  p <- substitute(p)

  pdf(paste0('results/', make.names(title_),".pdf"))
  res <- withVisible(eval(p))
  if (res$visible) print(res$value)
  if(sub_title) title(sub=title_)
  if(title) title(title_)
  dev.off()

  png(paste0('results/', make.names(title_),".png"),  width = 6, height = 6, units = 'in',res=600)
  res <- withVisible(eval(p))
  if (res$visible) print(res$value)
  if(sub_title) title(sub=title_)
  if(title) title(title_)
  dev.off()


  res <- withVisible(eval(p))
  if(sub_title) title(sub=title_)
  if(title) title(title_)
  if (res$visible) res$value else invisible(res$value)
}


stabilize_Anscombes <- function(count_matrix){
  # Use Anscombes approximation to variance stabilize Negative Binomial data
  # See https://f1000research.com/posters/4-1041 for motivation.
  # Assumes columns are samples, and rows are genes
  mu <- colMeans(count_matrix)
  vars <- matrixStats:::colVars(count_matrix)
  phi_hat <- optimize( function (phi) sum((mu + phi * mu^2-vars)^2), interval=c(1E-9,1E9))$minimum
  log2(count_matrix + 1. / (2 * phi_hat))
}

phi_hat <- function(mu, sigma_sq){
  # Use Anscombes approximation to variance stabilize Negative Binomial data
  # See https://f1000research.com/posters/4-1041 for motivation.
  # Assumes columns are samples, and rows are genes
  phi_hat <- optimize( function (phi) sum((mu + phi * mu^2-sigma_sq)^2), interval=c(1E-9,1E9))$minimum
  phi_hat
}

log_rcorr <- function(x, ...) {
  res <- Hmisc::rcorr(x, ...)
  h <- res$r
  npair <- res$n
  p <- as.integer(ncol(x))
  nam <- dimnames(x)[[2]]
  log_P <- matrix(log(2) + (
    pt(abs(h) * sqrt(npair - 2)/sqrt(1 - h * h), npair - 2, lower.tail = FALSE, log.p = TRUE)
  ), ncol = p)
  log_P[abs(h) == 1] <- 0
  diag(log_P) <- NA
  dimnames(log_P) <- list(nam, nam)
  res$log_P = log_P
  res
}

plot_corr_matrix <- function(X, ordering){
  lattice::levelplot(X[ordering,][,ordering], scales=list(draw=FALSE),
            xlab="gene", ylab="gene", col.regions = colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(201),
            at=seq(-1,1,0.01))
}

plot_corr_matrix <- function(X, ordering){
  lattice::levelplot(X[ordering,][,ordering], scales=list(draw=FALSE),
                     xlab="gene", ylab="gene", col.regions = viridis::viridis(201),
                     at=seq(-1,1,0.01))
}
cool_warm <- function(n) {
  colormap <- Rgnuplot:::GpdivergingColormap(seq(0,1,length.out=n),
                                             rgb1 = colorspace::sRGB( 0.230, 0.299, 0.754),
                                             rgb2 = colorspace::sRGB( 0.706, 0.016, 0.150),
                                             outColorspace = "sRGB")
  colormap[colormap>1] <- 1 # sometimes values are slightly larger than 1
  colormap <- grDevices::rgb(colormap[,1], colormap[,2], colormap[,3])
  colormap
}

plot_corr_matrix <- function(X, ordering){
  lattice::levelplot(X[ordering,][,ordering], scales=list(draw=FALSE),
                     xlab="gene", ylab="gene", col.regions = red_blue_diverging_colormap(500),
                     at=seq(-1,1,length.out=500))
}
