## Data generating processes.
dgp_NO <- function(n = 1000, probs = c(0.01, 0.1, 0.5, 0.9, 0.99),
  breaks = 20, het = TRUE, ...)
{
  ## Covariate.
  x <- runif(n, -3, 3)

  ## Parameters.
  mu <- sin(x)
  if(het)
    sigma <- exp(-1 + cos(x))
  else
    sigma <- 0.3

  ## Response
  y_cont <- rnorm(n, mean = mu, sd = sigma)

  ## Discretize the normal data into count categories.
  yr <- range(y_cont)
  breaks <- seq(-4, 4, length.out = breaks)
  y_count <- cut(y_cont, breaks = breaks, labels = FALSE, include.lowest = TRUE) - 1

  ## Combine.
  d <- data.frame("x" = x, "counts" = y_count, "num" = y_cont)

  ## Add quantiles.
  qu <- quc <- NULL
  for(j in probs) {
    qj <- qNO(j, mu = mu, sigma = sigma)
    qu <- cbind(qu, qj)
    quc <- cbind(quc, cut(qj, breaks = breaks, labels = FALSE, include.lowest = TRUE) - 1)
  }
  colnames(qu) <- colnames(quc) <- paste0(probs * 100, "%")
  d$quantiles <- qu
  d$count.quantiles <- quc

  return(d)
}


