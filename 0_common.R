
figure <- function(title, p){
  RESULT_FOLDER <- "results"
  ggsave(file.path(RESULT_FOLDER, paste0(make.names(title),".pdf")))
  p <- p + ggtitle(title)
  ggsave(file.path(RESULT_FOLDER, paste0(make.names(title),".png")))
  p
}

save_with_title <- function(title_){
  title(main = title_)
  dev.print(pdf, paste0('results/', make.names(title_), '.pdf'))
  dev.print(png, paste0('results/', make.names(title_), '.png'), width = 480, height = 480)
  invisible()
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
