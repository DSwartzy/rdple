ple0.jack.fn <- function(y, x, c, h, kernel="epanechnikov"){

  #getting the correct kernel function
  if (kernel=="epanechnikov"){mykernel="epanechnikov.fn"}
  else if (kernel=="triangular"){mykernel="triangular.fn"}
  else stop("You have chosen an invalid kernel.")
  #first getting the original estimate without anything deleted
  n <- length(x)
  l.kern.matrix <- kernelweight.fn(x.i=x, x.j=x, h=h, kernelname = mykernel, standardize.over = "x.j")
  f.vec <- as.numeric(x>c)

  f.fitted.ple0 <- t(l.kern.matrix)%*%f.vec
  f.resid.ple0 <- f.vec-f.fitted.ple0

  y.fitted.ple0 <- t(l.kern.matrix)%*%y
  y.resid.ple0 <- y-y.fitted.ple0

  tau.ple0 <- sum(f.resid.ple0^2)^(-1)*sum(f.resid.ple0*y.resid.ple0)

  df <- data.frame(f.ple0=f.resid.ple0, y.ple0=y.resid.ple0, index=seq(1,n))

  #getting Wu's jackknife with formulas that don't require a loop

  hat.ple0 <- (f.resid.ple0)%*%MASS::ginv(t(f.resid.ple0)%*%(f.resid.ple0))%*%t(f.resid.ple0)
  r.ple0 <- (diag(n) - hat.ple0)%*%(y.resid.ple0)

  se.ple0 <- sqrt(sum(df$f.ple0^2)^(-1)*sum(r.ple0^2/(1-diag(hat.ple0))*f.resid.ple0^2)*sum(df$f.ple0^2)^(-1))

  return(c(tau.ple0=tau.ple0, se.ple0=se.ple0))

}

ple1.jack.fn <- function(y, x, c, h, kernel="epanechnikov"){

  #getting the correct kernel function
  if (kernel=="epanechnikov"){mykernel="epanechnikov.fn"}
  else if (kernel=="triangular"){mykernel="triangular.fn"}
  else stop("You have chosen an invalid kernel.")

  #first getting the original estimate without anything deleted
  n <- length(x)
  l.ll.matrix <- sapply(1:n, LL.fn, x=x, h=h, kernelname=mykernel)
  f.vec <- as.numeric(x>c)

  f.fitted.ple1 <- t(l.ll.matrix)%*%f.vec
  f.resid.ple1 <- f.vec-f.fitted.ple1

  y.fitted.ple1 <- t(l.ll.matrix)%*%y
  y.resid.ple1 <- y-y.fitted.ple1

  tau.ple1 <- sum(f.resid.ple1^2)^(-1)*sum(f.resid.ple1*y.resid.ple1)

  df <- data.frame(f.ple1=f.resid.ple1, y.ple1=y.resid.ple1, index=seq(1,n))

  #getting Wu's jackknife with formulas that don't require a loop

  hat.ple1 <- (f.resid.ple1)%*%MASS::ginv(t(f.resid.ple1)%*%(f.resid.ple1))%*%t(f.resid.ple1)
  r.ple1 <- (diag(n) - hat.ple1)%*%(y.resid.ple1)

  se.ple1 <- sqrt(sum(df$f.ple1^2)^(-1)*sum(r.ple1^2/(1-diag(hat.ple1))*f.resid.ple1^2)*sum(df$f.ple1^2)^(-1))

  return(c(tau.ple1=tau.ple1, se.ple1=se.ple1))

}
