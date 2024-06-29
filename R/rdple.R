#' Title
#'
#' @param y A numeric vector.
#' @param x A numeric vector.
#' @param c A number.
#' @param h A number or a character string.
#' @param degree A number.
#' @param kernel A character string.
#' @param sparse.adj A character string.
#' @param sparse.r A number.
#'
#' @return A list.
#' @import stats
#'
rdple <- function(y, x, c, h="ska", degree=1, kernel="epanechnikov",
                     sparse.adj=FALSE, sparse.r=1){

  #checking that the inputs are vectors
  if(is.vector(x)==FALSE|is.vector(y)==FALSE) stop ("The x and y inputs must be vectors.")

  #checking for different lengths of x and y
  if (length(x)!=length(y)) stop ("The x and y inputs are of differing lengths.")

  #checking for missing values
  data <- data.frame(y,x)
  data2 <- na.omit(data)
  remove <- nrow(data)-nrow(data2)
  if (remove>0) warning(paste(remove, "observations with missing values have been removed."))

  x <- data2[,2]
  y <- data2[,1]

  #checking for c within the data
  n_below <- sum(as.numeric(x<c))
  n_above <- sum(as.numeric(x>=c))
  if (n_below==0|n_above==0) stop("Cutoff must be within the bounds of the running variable.")

  if (h=="ska"){h=ska_band(y=y, x=x, c=c, kernel=kernel)}
  else if (h<=0){stop("The bandwidth must be positive.")}
  else if (is.numeric(h)==FALSE) stop ("h must be a number or 'ska'.")

  #sparsity adjustment
  lgap <- sort(abs(x-c)[x<c])[1]
  rgap <- sort(abs(x-c)[x>=c])[1]
  mingap <- lgap + rgap

  if (sparse.adj==TRUE) {

    if (sparse.r<=0) stop ("The value of sparse.r must be positive.")
    else if (sparse.r%%1!=0)
      stop ("The value of sparse.r must be an integer.")
    else if (sparse.r>n_below|sparse.r>n_above)
      stop ("The value of sparse.r is larger than the number of observations above or below the cutoff.")

    #smallest nn to use as an epsilon
    neighbors <- dbscan::kNN(as.matrix(unique(x)), 1)
    epsilon <- min(neighbors$dist[,1])

    #gap on each side
    lgap <- sort(abs(x-c)[x<c])[sparse.r]
    rgap <- sort(abs(x-c)[x>=c])[sparse.r]
    mingap <- lgap + rgap

    #"minimum" bandwidth that works for our method
    if (h<mingap){h <- mingap + epsilon}
  }


  else if (sparse.adj==FALSE){
    if (mingap>=h) stop ("The chosen bandwidth does span the closest points on either side of the cutoff. Consider a sparsity adjustment.")
    if (sparse.r!=1) warning ("The chosen value of sparse.r was ignored since sparse.adj=FALSE.")
  }

  else stop ("The input for sparse.adj must be FALSE or TRUE.")

  if (degree==1){est=ple1.jack.fn(y=y, x=x, c=c, h=h, kernel=kernel)}
  else if (degree==0){est=ple0.jack.fn(y=y, x=x, c=c, h=h, kernel=kernel)}
  else stop("The degree must be 0 or 1.")
  return(list(tau=unname(est[1]), se=unname(est[2]), h=h))

}
