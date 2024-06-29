#' Title
#'
#' @param x A numeric vector.
#' @param c A number.
#'
#' @return A named numeric vector.
#' @import stats

diss <- function(x, c){

  #checking that the input is a vector
  if(is.vector(x)==FALSE) stop ("The x input must be a vector.")

  #checking for missing data
  data <- data.frame(x)
  data2 <- na.omit(data)
  remove <- nrow(data)-nrow(data2)
  if (remove>0) warning(paste(remove, "observations with missing values have been removed."))
  x <- data2$x

  #checking for c within the data
  n_below <- sum(as.numeric(x<c))
  n_above <- sum(as.numeric(x>=c))
  if (n_below==0|n_above==0) stop("Cutoff must be within the bounds of the running variable.")

  n <- length(x)
  s <- sd(x)
  scaled_iqr <- IQR(x)/1.34
  s_star <- min(s, scaled_iqr)
  h.rot <- 0.9*s_star*n^(-.2)
  m.l <- sum(as.numeric(((c-h.rot)<=x))&(x<c))
  m.r <- sum(as.numeric((x<(c+h.rot)))&(x>=c))
  m=m.l+m.r
  return(c(m=m, m.l=m.l, m.r=m.r, h.diss=h.rot))
}
