
figure <- function(title, p){
  RESULT_FOLDER <- "results"
  ggsave(file.path(RESULT_FOLDER, paste0(make.names(title),".pdf")))
  p <- p + ggtitle(title)
  ggsave(file.path(RESULT_FOLDER, paste0(make.names(title),".png")))
  p
}


figure1 <- function(title_, p, sub_title=FALSE){
  p <- substitute(p)

  pdf(paste0('results/', make.names(title_),".pdf"))
  eval(p)
  if(sub_title) title(sub=title_) else title(title_)
  dev.off()

  png(paste0('results/', make.names(title_),".png"))
  eval(p)
  if(sub_title) title(sub=title_) else title(title_)
  dev.off()

  ret <- eval(p)
  if(sub_title) title(sub=title_) else title(title_)
  invisible(ret)
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
