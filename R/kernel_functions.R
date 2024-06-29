#Epanechnikov function
epanechnikov.fn <- function(x) {
  0.75*(1-x^2)*(abs(x)<=1)
}

#Triangular function
triangular.fn <- function(x) {
  (1-abs(x))*(abs(x)<=1)
}

kernelweight.fn <- function(x.i, x.j, h, kernelname = "epanechnikov.fn", standardize.over="x.j") {
  #x.i is a vector of running variable values (or similar, like a single cutoff)
  #x.j is a vector of running variable values (or similar, like a single cutoff)
  #h is the scalar bandwidth
  #kernelname is the string name of the kernel function
  #standardize.over is an indicator of which vector the weights should be standardized over; options are:
  #"x.j" returns columns that sum to one
  #"x.i" returns rows that sum to one
  #"vector" returns a vector (not a matrix) that sums to one
  #"none" returns no normalization
  
  #this function takes two vectors and returns a **matrix** of normalized weights (even if x.i and x.j are scalars).
  #UNLESS the standardize.over option requests a vector, then it'll return a vector.
  
  #check to make sure some inputs are valid
  if(!(standardize.over) %in% c("x.j", "x.i", "vector", "none")) stop("you have asked for an invalid standardization.")
  if(standardize.over == "vector" & length(x.i)>1 & length(x.j)>1) stop("the vector option is only valid if at least one inputs are scalar.")
  if((standardize.over =="x.j" & length(x.i)==1) | (standardize.over =="x.i" & length(x.j)==1)) warning("I think you have chosen the incorrect standardization margin.") 
  
  #first, calculate the distance between each pair of elements in the two vectors.
  mydistance <- outer(x.i/h, x.j/h, "-") #this is a matrix where rows correspond to x.i an columns correspond to x.j
  #here we calculate unscaled weights using the generic name for the kernel function -- thus we can swap it out easily later.
  myunscaledweight <- do.call(kernelname, list(x=mydistance)) #this is a matrix of the same dimensions as mydistance
  
  #now we need to standardize
  if(standardize.over == "x.j")  myweights <- prop.table(myunscaledweight, margin = 2)
  if(standardize.over == "x.i")  myweights <- prop.table(myunscaledweight, margin = 1)
  if(standardize.over == "vector") myweights <- as.vector(myunscaledweight)/sum(as.vector(myunscaledweight))
  if(standardize.over == "none")  myweights <- myunscaledweight
  
  return(myweights)
}

#Weight function for local linear regression
LL.fn <- function(i, x, h, n, kernelname= "epanechnikov.fn"){
  #requires function kernelweight.fn
  
  if(!("MASS" %in% .packages(all.available=TRUE))) stop("You must install the MASS package.")
  
  w.xi <- as.vector(kernelweight.fn(x[i], x, h, standardize.over="none", kernelname = kernelname))
  
  #create the x design matrix
  x.xi <- cbind(rep(1,n), x-x[i])
  
  #create the e vector
  e.vec <- c(1,0)
  
  #get the results
  return((w.xi*x.xi)%*%MASS::ginv(t(x.xi)%*%(w.xi*x.xi))%*%e.vec)
}
