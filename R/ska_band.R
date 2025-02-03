#' Swartzentruber-KAizar bandwidth algorithm
#'
#' Calculates the Swartzentruber-KAizar (SKA) bandwidth value for an RD design.
#'
#' @param y the response variable.
#' @param x the running variable.
#' @param c the cutoff.
#' @param kernel specifies the kernel function used for the local polynomial weights. Options are `epanechnikov` (default) or `triangular`.
#' @return The SKA bandwidth.
#' @export
#' @import stats
#'

ska_band <- function(y, x, c, kernel="epanechnikov") {

  #checking that the inputs are vectors
  if(is.vector(x)==FALSE|is.vector(y)==FALSE) stop ("The x and y inputs must be vectors.")

  #checking for different lengths of x and y
  if (length(x)!=length(y)) stop ("The x and y inputs are of differing lengths.")

  #checking for missing data
  data <- data.frame(y,x)
  data2 <- na.omit(data)
  remove <- nrow(data)-nrow(data2)
  if (remove>0) warning(paste(remove, "observations with missing values have been removed."))

  x <- data2$x
  y <- data2$y

  #checking for c within the data
  n_below <- sum(as.numeric(x<c))
  n_above <- sum(as.numeric(x>=c))
  if (n_below==0|n_above==0) stop("Cutoff must be within the bounds of the running variable.")
  else if (n_below<3|n_above<3)
    stop("At least 3 values on both sides of the cutoff are necessary to calculate the bandwidth.")

  if(kernel=="epanechnikov"){kernelconst=epaconst}
  else if (kernel=="triangular"){kernelconst=triconst}
  else stop("You have chosen an invalid kernel.")

  v.k <- kernelconst[29]
  a.p2 <- kernelconst[1:3]
  b.p2 <- kernelconst[13:15]
  a.p3 <- kernelconst[4:7]
  b.p3 <- kernelconst[16:19]
  a.p4 <- kernelconst[8:12]
  b.p4 <- kernelconst[20:24]
  b.b.1 <- kernelconst[25]
  b.b.2 <- kernelconst[26]
  b.b.3 <- kernelconst[27]
  b.b.4 <- kernelconst[28]

  n <- length(x)

  #using options from the ks package

  h0.pilot1 <- ks::hpi(x, deriv.order=0)
  h1.pilot1 <- ks::hpi(x, deriv.order=1)
  h2.pilot1 <- ks::hpi(x, deriv.order=2)

  f0 <- ks::kdde(x, h=h0.pilot1, deriv.order=0, eval.points=c(c))$estimate
  f1 <- ks::kdde(x, h=h1.pilot1, deriv.order=1, eval.points=c(c))$estimate
  f2 <- ks::kdde(x, h=h2.pilot1, deriv.order=2, eval.points=c(c))$estimate

  #getting the variance estimates based on the nearest neighbor idea from CCT, with 3 neighbors

  x.u <- unique(x)
  y.u <- y[duplicated(x)==FALSE]
  x.right <- x.u[x.u>=c]
  y.right <- y.u[x.u>=c]
  right.index <- order(x.right)[1:3]
  x.right.nn <- x.right[right.index]
  y.right.nn <- y.right[right.index]
  x.left <- x[x<c]
  y.left <- y[x<c]
  left.index <- order(x.left, decreasing=TRUE)[1:3]
  x.left.nn <- x.left[left.index]
  y.left.nn <- y.left[left.index]
  sigma2.hat.right <- sum((y.right.nn-mean(y.right.nn))^2)/(length(y.right.nn)-1)
  sigma2.hat.left <- sum((y.left.nn-mean(y.left.nn))^2)/(length(y.left.nn)-1)

  #get the mean function derivative estimates

  #using a global cubic function to estimate the third derivative of the mean and sigma (here we assume equal variance)
  pilot.global3.lm <- lm(y ~ I(as.numeric((x-c)>0)) + I(x-c) + I((x-c)^2) + I((x-c)^3))
  m3.hat.1 <- 6*pilot.global3.lm$coef[[5]]
  sigma.pilot3.hat <- summary(pilot.global3.lm)$sigma

  #using a global quartic function to estimate the fourth derivative of the mean and sigma (here we assume equal variance)
  pilot.global4.lm <- lm(y ~ I(as.numeric((x-c)>0)) + I(x-c) + I((x-c)^2) + I((x-c)^3) + I((x-c)^4))
  m4.hat <- 24*pilot.global4.lm$coef[[6]]
  sigma.pilot4.hat <- summary(pilot.global4.lm)$sigma

  #using a global quintic function to estimate the fifth derivative of the mean and sigma (here we assume equal variance)
  pilot.global5.lm <- lm(y ~ I(as.numeric((x-c)>0)) + I(x-c) + I((x-c)^2) + I((x-c)^3) + I((x-c)^4) + I((x-c)^5))
  m5.hat <- 120*pilot.global5.lm$coef[[7]]
  sigma.pilot5.hat <- summary(pilot.global5.lm)$sigma

  #plugging in our estimates to the pilot bandwidth formula
  #when v=1 I can use the cubic (p=2)
  v <- 1
  h.pilot.v1 <- (((2*v+1)*(factorial(2+1))^2*(a.p2[v+1])*sigma.pilot3.hat^2)/(2*(3-v)*(b.p2[v+1]^2)*m3.hat.1^2*n*f0))^(1/7)

  #when v=2 I can use the quartic (p=3)
  v <- 2
  h.pilot.v2 <- (((2*v+1)*(factorial(3+1))^2*(a.p3[v+1])*sigma.pilot4.hat^2)/(2*(4-v)*(b.p3[v+1]^2)*m4.hat^2*n*f0))^(1/9)

  #when v=3 I can use the quintic (p=4)
  v <- 3
  h.pilot.v3 <- (((2*v+1)*(factorial(4+1))^2*(a.p4[v+1])*sigma.pilot5.hat^2)/(2*(5-v)*(b.p4[v+1]^2)*m5.hat^2*n*f0))^(1/11)


  #keeping only values within these bandwidths of the cutoff
  x.keep.v1 <- x[x>(c-h.pilot.v1) & x<(c+h.pilot.v1)]
  x.keep.v2 <- x[x>(c-h.pilot.v2) & x<(c+h.pilot.v2)]
  x.keep.v3 <- x[x>(c-h.pilot.v3) & x<(c+h.pilot.v3)]
  y.keep.v1 <- y[x>(c-h.pilot.v1) & x<(c+h.pilot.v1)]
  y.keep.v2 <- y[x>(c-h.pilot.v2) & x<(c+h.pilot.v2)]
  y.keep.v3 <- y[x>(c-h.pilot.v3) & x<(c+h.pilot.v3)]

  #running a v+1 degree polynomial jump for v=1,2,3
  local.v1.lm <- lm(y.keep.v1 ~ I(as.numeric((x.keep.v1-c)>0)) + I(x.keep.v1-c) + I((x.keep.v1-c)^2))
  local.v2.lm <- lm(y.keep.v2 ~ I(as.numeric((x.keep.v2-c)>0)) + I(x.keep.v2-c) + I((x.keep.v2-c)^2) + I((x.keep.v2-c)^3))
  local.v3.lm <- lm(y.keep.v3 ~ I(as.numeric((x.keep.v3-c)>0)) + I(x.keep.v3-c) + I((x.keep.v3-c)^2) + I((x.keep.v3-c)^3) + I((x.keep.v3-c)^4))

  m1.hat <- local.v1.lm$coefficient[[3]]
  m2.hat <- 2*local.v2.lm$coefficient[[4]]
  m3.hat <- 6*local.v3.lm$coefficient[[5]]

  #calculating b.b (will change with a different kernel)
  g2 <- m1.hat*f1+m2.hat*f0/2
  g2.prime <- m1.hat*f2+m2.hat*f1+(m2.hat*f1+m3.hat*f0)/2
  b.b <- 2*b.b.1*(f0*b.b.2)^(-1)*(f1/f0*g2*b.b.3-g2.prime*b.b.4)

  #calculating the final bandwidth
  h.ska <- (v.k*(sigma2.hat.left+sigma2.hat.right)/n/b.b^2/24/f0)^(1/7)


  return(h.ska)

}
