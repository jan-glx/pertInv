#' @export
deconv <-function (tau, X, y, Q, P, n = 40, family = c("Poisson", "PoissonExp","PoissonExpLog", "Normal", "Binomial"), ignoreZero = TRUE, deltaAt = NULL, c0 = 1, scale = TRUE,
                    pDegree = 5, aStart = 1,  ...)
{
  family <- match.arg(family)
  if (missing(Q) && missing(P)) {
    if (family == "Poisson") {
      m <- length(tau)
      if (ignoreZero) {
        supportOfX <- seq_len(n)
        P <- sapply(tau, function(lam) dpois(x = supportOfX,
                                             lambda = lam)/(1 - exp(-lam)))
      }
      else {
        supportOfX <- seq.int(from = 0, to = n - 1)
        P <- sapply(tau, function(w) dpois(x = supportOfX,
                                           lambda = w))
      }
      if (missing(y)) {
        y <- sapply(supportOfX, function(i) sum(X ==
                                                  i))
      }
      Q <- cbind(1, scale(splines::ns(tau, pDegree), center = TRUE,
                          scale = FALSE))
      Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w)))
    } else if (family == "PoissonExp") {

      if (missing(tau)) {
        tau <- seq(0, max(y)/mean(X), length.out = n+1)[-1]
        m <- length(tau)
      }
      if (ignoreZero) {
        P <- sapply(tau, function(lam) dpois(x = y, lambda = lam*X)/(1 - exp(-lam*X)))
      } else {
        P <- sapply(tau, function(lam) dpois(x = y, lambda = lam*X))
      }

      Q <- cbind(1, scale(splines::ns(tau, pDegree), center = TRUE, scale = FALSE))
      Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w)))
      #browser()
      y <- 1
    } else if (family == "PoissonExpLog") {

      if (missing(tau)) {
        tau <- seq(log(min(y+1)/10), log(max(y)/mean(X)), length.out = n)
        m <- length(tau)
      }
      if (ignoreZero) {
        P <- sapply(tau, function(lam) dpois(x = y, lambda = exp(lam)*X)/(1 - exp(-exp(lam)*X)))
      } else {
        P <- sapply(tau, function(lam) dpois(x = y, lambda = exp(lam)*X))
      }

      Q <- cbind(1, scale(splines::ns(tau, pDegree), center = TRUE, scale = FALSE))
      Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w)))
      #browser()
      y <- 1
    }
    else if (family == "Normal") {
      m <- length(tau)
      r <- round(range(X), digits = 1)
      xBin <- seq(from = r[1], to = r[2], length.out = n)
      xBinDropFirst <- xBin[-1]
      xBinDropLast <- xBin[-length(xBin)]
      P <- sapply(tau, function(x) pnorm(q = xBinDropFirst,
                                         mean = x) - pnorm(q = xBinDropLast, mean = x))
      intervals <- findInterval(X, vec = xBin)
      y <- sapply(seq_len(n - 1), function(w) sum(intervals ==
                                                    w))
      if (scale) {
        Q1 <- scale(splines::ns(tau, pDegree), center = TRUE,
                    scale = FALSE)
        Q1 <- apply(Q1, 2, function(w) w/sqrt(sum(w *
                                                    w)))
      }
      if (!is.null(deltaAt)) {
        I0 <- as.numeric(abs(tau - deltaAt) < 1e-10)
        Q <- cbind(I0, Q1)
      }
      else {
        Q <- Q1
      }
    }
    else {
      m <- length(tau)
      Q <- splines::ns(tau, pDegree)
      if (scale) {
        Q <- scale(Q, center = TRUE, scale = FALSE)
        Q <- apply(Q, 2, function(w) w/sqrt(sum(w * w)))
      }
      P <- sapply(tau, function(w) dbinom(X[, 2], size = X[,
                                                           1], prob = w))
      y <- 1
    }
  }
  else {
    if (!missing(X) || missing(y) || missing(P) || missing(Q)) {
      stop("P, Q, and y (but not X) must be specified together!")
    }
  }
  p <- ncol(Q)
  pGiven <- length(aStart)
  if (pGiven == 1) {
    aStart <- rep(aStart, p)
  }
  else {
    if (pGiven != p)
      stop(sprintf("Wrong length (%d) for initial parameter, expecting length (%d)",
                   pGiven, p))
  }
  statsFunction <- function(a) {
    g <- as.vector(exp(Q %*% a))
    g <- g/sum(g)
    G <- cumsum(g)
    f <- as.vector(P %*% g)
    yHat <- sum(y) * f
    Pt <- P/f
    W <- g * (t(Pt) - 1)
    qw <- t(Q) %*% W
    ywq <- (yHat * t(W)) %*% Q
    I1 <- qw %*% ywq
    aa <- sqrt(sum(a^2))
    sDot <- c0 * a/aa
    sDotDot <- (c0/aa) * (diag(length(a)) - outer(a, a)/aa^2)
    R <- sum(diag(sDotDot))/sum(diag(I1))
    I2 <- solve(I1 + sDotDot)
    bias <- as.vector(-I2 %*% sDot)
    Cov <- I2 %*% (I1 %*% t(I2))
    Dq <- (diag(g) - outer(g, g)) %*% Q
    bias.g <- Dq %*% bias
    Cov.g <- Dq %*% Cov %*% t(Dq)
    se.g <- diag(Cov.g)^0.5
    D <- diag(length(tau))
    D[lower.tri(D)] <- 1
    Cov.G <- D %*% (Cov.g %*% t(D))
    se.G <- diag(Cov.G)^0.5
    mat <- cbind(tau, g, se.g, G, se.G, bias.g)
    colnames(mat) = c("theta", "g", "SE.g", "G", "SE.G",
                      "Bias.g")
    list(S = R, cov = Cov, cov.g = Cov.g, mat = mat)
  }
  loglik <- function(a) {
    g <- exp(Q %*% a)
    g <- as.vector(g/sum(g))
    f <- as.vector(P %*% g)
    value <- -sum(y * log(f)) + c0 * sum(a^2)^0.5
    Pt <- P/f
    W <- g * (t(Pt) - 1)
    qw <- t(Q) %*% W
    aa <- sqrt(sum(a^2))
    sDot <- c0 * a/aa
    attr(value, "gradient") <- if (identical(y, 1)) -rowSums(qw) + sDot else  -(qw %*% y) + sDot
    value
  }
  result <- stats::nlm(f = loglik, p = aStart, gradtol = 1e-10,
                       ...)
  if (result$code>2) warning("nlm termination code ", result$code)
  mle <- result$estimate
  stats <- statsFunction(mle)
  list(mle = mle, Q = Q, P = P, S = stats$S, cov = stats$cov,
       cov.g = stats$cov.g, stats = stats$mat, loglik = loglik,
       statsFunction = statsFunction)
}
