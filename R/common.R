
#' @export
figure <- function(title, p, sub_title = if(exists("data_set")) data_set else character(0), ...){
  RESULT_FOLDER <- "results"
  file_name <- paste0(make.names(title), "-", make.names(sub_title))
  cowplot::ggsave(file.path(RESULT_FOLDER, paste0(file_name,".pdf")), p, ...)#, device=cairo_pdf)
  cowplot::ggsave(file.path(RESULT_FOLDER, paste0(file_name,".png")), p, ..., type = "cairo")
  p <- p + ggplot2::ggtitle(title, sub_title)
  cowplot::ggsave(file.path(RESULT_FOLDER, paste0(file_name,"_titled.png")), p, ..., type = "cairo")
  p
}

#' Saves last figure as png and as pdf
#'
#'
#'
#' @param title_ filename for save optional title for the plot
#' @param p unevaluated expression that performs the plot
#' @param sub_title \code{title_} will be added as sub title to the plot if set
#' @param title \code{title_} will be added as title to the plot if set
#' @param width numeric, width of the saved plots
#' @param height numeric, height of the saved plots
#' @param units character, unit of \code{width} and \code{height} values
#' @param res numeric, resolution of the saved plots
#' @param ... further arguments passed to \code{png()} and \code{pdf()}
#'
#' @return \code{p} evaluated
#'
#'
#' @importFrom grDevices dev.off pdf png
#' @export
figure1 <- function(title_, p, sub_title=FALSE, title=!sub_title, width=4, height=4, units = 'in', res=600, ...){
  p <- substitute(p)

  pdf(paste0('results/', make.names(title_),".pdf"), width=width, height=height, ...)
  ret <- withVisible(eval(p))
  if (ret$visible) print(ret$value)
  #if(sub_title) title(sub=title_)
  #if(title) title(title_)
  dev.off()

  png(paste0('results/', make.names(title_),".png"), width=width, height=height, units=units, res=res, type="cairo", ...)
  ret <- withVisible(eval(p))
  if (ret$visible) print(ret$value)
  if(sub_title) title(sub=title_)
  if(title) title(title_)
  dev.off()


  ret <- withVisible(eval(p))
  if(sub_title) title(sub=title_)
  if(title) title(title_)
  if (ret$visible) ret$value else invisible(ret$value)
}

#' @export
stabilize_Anscombes <- function(count_matrix){
  # Use Anscombes approximation to variance stabilize Negative Binomial data
  # See https://f1000research.com/posters/4-1041 for motivation.
  # Assumes columns are samples, and rows are genes
  mu <- colMeans(count_matrix)
  sigma_sq <- matrixStats::colVars(count_matrix)
  phi_hat <- phi_hat(mu, sigma_sq)
  log2(count_matrix + 1. / (2 * phi_hat))
}

#' @export
phi_hat <- function(mu, sigma_sq){
  # Use Anscombes approximation to variance stabilize Negative Binomial data
  # See https://f1000research.com/posters/4-1041 for motivation.
  # Assumes columns are samples, and rows are genes
  phi_hat <- stats::optimize( function (phi) sum((mu + phi * mu^2-sigma_sq)^2), interval=c(1E-9,1E9))$minimum
  phi_hat
}

#' @export
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

#' @export
plot_corr_matrix <- function(X, ordering=seq_len(dim(X)[1])){
  lattice::levelplot(X[ordering,][,ordering], scales=list(draw=FALSE),
            xlab="gene", ylab="gene", col.regions = colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(201),
            at=seq(-1,1,0.01))
}

#' @export
plot_matrix <- function(X){
  lattice::levelplot(X, scales=list(draw=FALSE),
                     xlab="gene", ylab="gene", col.regions = viridis::viridis(101),
                     at=seq(0,1,0.01))
}

#  http://jfly.iam.u-tokyo.ac.jp/color/:
#' @export
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
#' @export
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#' @export
cool_warm <- function(n) {
  colormap <- Rgnuplot::GpdivergingColormap(seq(0,1,length.out=n),
                                             rgb1 = colorspace::sRGB( 0.230, 0.299, 0.754),
                                             rgb2 = colorspace::sRGB( 0.706, 0.016, 0.150),
                                             outColorspace = "sRGB")
  colormap[colormap>1] <- 1 # sometimes values are slightly larger than 1
  colormap <- grDevices::rgb(colormap[,1], colormap[,2], colormap[,3])
  colormap
}

#' @export
plot_corr_matrix <- function(X, ordering){
  lattice::levelplot(X[ordering,][,ordering], scales=list(draw=FALSE),
                     xlab="gene", ylab="gene", col.regions = cool_warm(500),
                     at=seq(-1,1,length.out=500))
}

#' @export
is_dev <- function(){
  exists("development_mode") && development_mode
}


GeomSplitViolin <- ggplot2::ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                              1))
    quantiles <- create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

#' @export
geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}


#' @export
cross_entropy_loss <-  function(p, label) -sum(log((!label)+(2*label-1)*p))/length(label)

#' @export
cross_entropy_loss_llr = function(LLR, label) {
  mean(log_sum_exp(0, LLR*(1-2*label)))
}


#' @export
mean_sd <- function (x, mult = 1)
{
  std <- mult * stats::sd(x, na.rm=TRUE)
  mean <- mean(x, na.rm=TRUE)
  data.frame(y = mean, ymin = mean - std, ymax = mean + std)
}

#' @export
cell_loglikelihoods <- function(res, sigma) apply(res, MARGIN=1, function(x) sum(dnorm(x, sd=sigma,log = TRUE)))

#' @export
quantile_ci <-  function(y, p=0.5, conf.level = 0.95, alpha = 1-conf.level) {
  rbindlist(lapply(p, function(p) {
    m <-  quantile(y, probs=p)
    n <-  length(y)
    lower_i <- qbinom(alpha/2, n, p)
    upper_i <- 1+qbinom(alpha/2, n, p, lower.tail=FALSE)
    ymin <- -Inf
    ymax <- Inf
    if (lower_i>0 && upper_i<n) {
      qs <- sort(y, partial=c(lower_i,upper_i))[c(lower_i,upper_i)]
      ymin <-  qs[1]
      ymax <-  qs[2]
    } else if (lower_i>0) {
      ymin <- sort(y, partial=lower_i)[lower_i]
    } else if (upper_i<n){
      ymax <- sort(y, partial=upper_i)[upper_i]
    }
    data.table(ymin, y=m, ymax, p = p)
  }))
}

#' @export
PositionJitterNormal <- ggplot2::ggproto("PositionJitterNormal", PositionJitter, compute_layer = function(self, data, params, panel){
  trans_x <- if (params$width > 0)
    function(x) x+rnorm(length(x), sd= params$width)
  trans_y <- if (params$height > 0)
    function(x) x+rnorm(length(x), sd= params$height)
  transform_position(data, trans_x, trans_y)
})

#' @export
position_jitter_normal <- function (width = NULL, height = NULL)
{
  ggplot2::ggproto(NULL, PositionJitterNormal, width = width, height = height)
}


#' @export
geom_jitter_normal <-function (mapping = NULL, data = NULL, stat = "identity", position = "jitter",
                               ..., width = NULL, height = NULL, na.rm = FALSE, show.legend = NA,
                               inherit.aes = TRUE)
{
  if (!missing(width) || !missing(height)) {
    if (!missing(position)) {
      stop("Specify either `position` or `width`/`height`",
           call. = FALSE)
    }
    position <- position_jitter_normal(width = width, height = height)
  }
  ggplot2::layer(data = data, mapping = mapping, stat = stat, geom = GeomPoint,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ...))
}

#' @export
sample.dt <- function(x, size, replace= FALSE, prob = NULL) {
  x[sample(nrow(x), size, replace, prob)]
}


#' @export
index_first_dim <- function(x, i) {
  sx <- quote(x[i])
  for (j in seq_len(length(dim(x))-1)) {
    sx[[j+3]] <- quote(x[,])[[3]] #add empty symbol
  }
  eval(sx)
}

#' @export
index_sample <-  function(samples, i) {
  lapply(samples, index_first_dim, i=i)
}

#' @export
index_samples <-  function(samples, is) {
  lapply(is, function(i) index_sample(samples, i))
}

#' @export
log_sum_exp <- function(x, y=NULL) {
  if (missing(y)){
    xm <- max(x)
    log(sum(exp(x - xm))) + xm
  } else {
    xm <- pmax(x, y)
    res <- log(rowSums(cbind(as.vector(exp(x - xm)), as.vector(exp(y - xm))))) + as.vector(xm)
    d <- pmax(dim(x), dim(y))
    if (length(d)) dim(res) <- d
    res
  }
}

#' @export
log_mean_exp <- function(x) {
  log_sum_exp(x)-log(length(x))
}

#' @export
#' @author gvrocha
#' @source https://stackoverflow.com/a/33179099/1870254
log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)

    minor_breaks =log10(seq(1, 9, by=1))
    breaks = rep(major_breaks, each=length(minor_breaks)) +
             rep(minor_breaks, times = n_major)
    if(sum( (breaks < max(log10(x))) & (breaks > min(log10(x)))) >5 ){
      breaks <-  breaks[rep(c(T,T,T,F,T,F,T,F,F),n_major)]
    }
    return(10^(breaks))
  }
}

#' @export
binom_mean_ci_jeffreys <-  function(x, ...) as.data.table(binom::binom.bayes(sum(x),length(x), ...))[,.(ymin=lower,ymax=upper,y=mean)]

