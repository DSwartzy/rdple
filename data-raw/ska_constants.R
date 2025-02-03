#This script generates and saves constants for the skaband function
#Currently uses Epanechnikov and triangular kernels

#Epanechnikov kernel====
#calculating constants for the epanechnikov kernel
s.0 <- integrate(function(x){x^0*3/4*(1-x^2)}, -1, 1)$value
s.1 <- integrate(function(x){x^1*3/4*(1-x^2)}, -1, 1)$value
s.2 <- integrate(function(x){x^2*3/4*(1-x^2)}, -1, 1)$value
s.3 <- integrate(function(x){x^3*3/4*(1-x^2)}, -1, 1)$value
s.4 <- integrate(function(x){x^4*3/4*(1-x^2)}, -1, 1)$value
s.5 <- integrate(function(x){x^5*3/4*(1-x^2)}, -1, 1)$value
s.6 <- integrate(function(x){x^6*3/4*(1-x^2)}, -1, 1)$value
s.7 <- integrate(function(x){x^7*3/4*(1-x^2)}, -1, 1)$value
s.8 <- integrate(function(x){x^8*3/4*(1-x^2)}, -1, 1)$value
s.9 <- integrate(function(x){x^9*3/4*(1-x^2)}, -1, 1)$value
v.0 <- integrate(function(x){x^0*(3/4*(1-x^2))^2}, -1, 1)$value
v.1 <- integrate(function(x){x^1*(3/4*(1-x^2))^2}, -1, 1)$value
v.2 <- integrate(function(x){x^2*(3/4*(1-x^2))^2}, -1, 1)$value
v.3 <- integrate(function(x){x^3*(3/4*(1-x^2))^2}, -1, 1)$value
v.4 <- integrate(function(x){x^4*(3/4*(1-x^2))^2}, -1, 1)$value
v.5 <- integrate(function(x){x^5*(3/4*(1-x^2))^2}, -1, 1)$value
v.6 <- integrate(function(x){x^6*(3/4*(1-x^2))^2}, -1, 1)$value
v.7 <- integrate(function(x){x^7*(3/4*(1-x^2))^2}, -1, 1)$value
v.8 <- integrate(function(x){x^8*(3/4*(1-x^2))^2}, -1, 1)$value

#when p=2
s <- matrix(c(s.0, s.1, s.2, s.1, s.2, s.3, s.2, s.3, s.4), nrow=3)
s.star <- matrix(c(v.0, v.1, v.2, v.1, v.2, v.3, v.2, v.3, v.4), nrow=3)

a.p2 <- diag(solve(s)%*%s.star%*%solve(s))
b.p2 <- as.vector(solve(s)%*%c(s.3, s.4, s.5))

#when p=3
s <- matrix(c(s.0, s.1, s.2, s.3, s.1, s.2, s.3, s.4, 
              s.2, s.3, s.4, s.5, s.3, s.4, s.5, s.6), nrow=4)
s.star <- matrix(c(v.0, v.1, v.2, v.3, v.1, v.2, v.3, v.4, 
                   v.2, v.3, v.4, v.5, v.3, v.4, v.5, v.6), nrow=4)

a.p3 <- diag(solve(s)%*%s.star%*%solve(s))
b.p3 <- as.vector(solve(s)%*%c(s.4, s.5, s.6, s.7))

#when p=4
s <- matrix(c(s.0, s.1, s.2, s.3, s.4, s.1, s.2, s.3, s.4, s.5,
              s.2, s.3, s.4, s.5, s.6, s.3, s.4, s.5, s.6, s.7,
              s.4, s.5, s.6, s.7, s.8), nrow=5)
s.star <- matrix(c(v.0, v.1, v.2, v.3, v.4, v.1, v.2, v.3, v.4, v.5,
                   v.2, v.3, v.4, v.5, v.6, v.3, v.4, v.5, v.6, v.7,
                   v.4, v.5, v.6, v.7, v.8), nrow=5)

a.p4 <- diag(solve(s)%*%s.star%*%solve(s))
b.p4 <- as.vector(solve(s)%*%c(s.5, s.6, s.7, s.8, s.9))

#now getting v.k which is called Cp1 in the dissertation
K_0.fn <- function(x){.25*(x^3-3*x+2)}
K_1.fn <- function(x){(3/16)*(x^4-2*x^2+1)}
K_2.fn <- function(x){(1/20)*(3*x^5-5*x^3+2)}
L.fn <- function(x){(1/320)*(-45*x^2-24*x+40)}
square <- function(x){x^2}

v.k <- integrate(function(x){square(K_0.fn(x))}, 0,1)$value^(-2)*
integrate(function(x){square(K_0.fn(x)+L.fn(x)-L.fn(-x))}, 0,1)$value

#constants for the b.b expression
b.b.1 <- K_2.fn(0)
b.b.2 <- integrate(function(x){square(K_0.fn(x))}, 0,1)$value
b.b.3 <- integrate(function(x){K_1.fn(x)}, 0,1)$value
b.b.4 <- integrate(function(x){K_0.fn(x)*x}, 0,1)$value

#saving as a vector
kernelconst <- c(a.p2, a.p3, a.p4, b.p2, b.p3, b.p4, b.b.1, b.b.2, b.b.3, b.b.4, v.k)
save(kernelconst, file="epaconst.RData")

#Triangular kernel====
s.0 <- integrate(function(x){x^0*(1-abs(x))}, -1, 1)$value
s.1 <- integrate(function(x){x^1*(1-abs(x))}, -1, 1)$value
s.2 <- integrate(function(x){x^2*(1-abs(x))}, -1, 1)$value
s.3 <- integrate(function(x){x^3*(1-abs(x))}, -1, 1)$value
s.4 <- integrate(function(x){x^4*(1-abs(x))}, -1, 1)$value
s.5 <- integrate(function(x){x^5*(1-abs(x))}, -1, 1)$value
s.6 <- integrate(function(x){x^6*(1-abs(x))}, -1, 1)$value
s.7 <- integrate(function(x){x^7*(1-abs(x))}, -1, 1)$value
s.8 <- integrate(function(x){x^8*(1-abs(x))}, -1, 1)$value
s.9 <- integrate(function(x){x^9*(1-abs(x))}, -1, 1)$value
v.0 <- integrate(function(x){x^0*((1-abs(x)))^2}, -1, 1)$value
v.1 <- integrate(function(x){x^1*((1-abs(x)))^2}, -1, 1)$value
v.2 <- integrate(function(x){x^2*((1-abs(x)))^2}, -1, 1)$value
v.3 <- integrate(function(x){x^3*((1-abs(x)))^2}, -1, 1)$value
v.4 <- integrate(function(x){x^4*((1-abs(x)))^2}, -1, 1)$value
v.5 <- integrate(function(x){x^5*((1-abs(x)))^2}, -1, 1)$value
v.6 <- integrate(function(x){x^6*((1-abs(x)))^2}, -1, 1)$value
v.7 <- integrate(function(x){x^7*((1-abs(x)))^2}, -1, 1)$value
v.8 <- integrate(function(x){x^8*((1-abs(x)))^2}, -1, 1)$value

#when p=2
s <- matrix(c(s.0, s.1, s.2, s.1, s.2, s.3, s.2, s.3, s.4), nrow=3)
s.star <- matrix(c(v.0, v.1, v.2, v.1, v.2, v.3, v.2, v.3, v.4), nrow=3)

a.p2 <- diag(solve(s)%*%s.star%*%solve(s))
b.p2 <- as.vector(solve(s)%*%c(s.3, s.4, s.5))

#when p=3
s <- matrix(c(s.0, s.1, s.2, s.3, s.1, s.2, s.3, s.4, 
              s.2, s.3, s.4, s.5, s.3, s.4, s.5, s.6), nrow=4)
s.star <- matrix(c(v.0, v.1, v.2, v.3, v.1, v.2, v.3, v.4, 
                   v.2, v.3, v.4, v.5, v.3, v.4, v.5, v.6), nrow=4)

a.p3 <- diag(solve(s)%*%s.star%*%solve(s))
b.p3 <- as.vector(solve(s)%*%c(s.4, s.5, s.6, s.7))

#when p=4
s <- matrix(c(s.0, s.1, s.2, s.3, s.4, s.1, s.2, s.3, s.4, s.5,
              s.2, s.3, s.4, s.5, s.6, s.3, s.4, s.5, s.6, s.7,
              s.4, s.5, s.6, s.7, s.8), nrow=5)
s.star <- matrix(c(v.0, v.1, v.2, v.3, v.4, v.1, v.2, v.3, v.4, v.5,
                   v.2, v.3, v.4, v.5, v.6, v.3, v.4, v.5, v.6, v.7,
                   v.4, v.5, v.6, v.7, v.8), nrow=5)

a.p4 <- diag(solve(s)%*%s.star%*%solve(s))
b.p4 <- as.vector(solve(s)%*%c(s.5, s.6, s.7, s.8, s.9))


#now getting v.k which is called Cp1 in the dissertation
K_0.fn <- function(x){.5*(x-1)^2}
K_1.fn <- function(x){(1/6)*(x-1)^2*(2*x+1)}
K_2.fn <- function(x){(1/12)*(3*x^4-4*x^3+1)}
L.fn <- function(x){.125-(1/6)*x}
square <- function(x){x^2}

v.k <- integrate(function(x){square(K_0.fn(x))}, 0,1)$value^(-2)*
  integrate(function(x){square(K_0.fn(x)+L.fn(x)-L.fn(-x))}, 0,1)$value

#constants for the b.b expression
b.b.1 <- K_2.fn(0)
b.b.2 <- integrate(function(x){square(K_0.fn(x))}, 0,1)$value
b.b.3 <- integrate(function(x){K_1.fn(x)}, 0,1)$value
b.b.4 <- integrate(function(x){K_0.fn(x)*x}, 0,1)$value

#saving as a vector
kernelconst <- c(a.p2, a.p3, a.p4, b.p2, b.p3, b.p4, b.b.1, b.b.2, b.b.3, b.b.4, v.k)
save(kernelconst, file="triconst.RData")
