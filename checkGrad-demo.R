##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
library(R.utils)
sourceDirectory("lib", modifiedOnly=FALSE) # Source the 'lib' directory

f <- c(4, .5, .7, .2, .9)
m <- c(10, 11, 10, 11, 10)
k <- c(5, 7, 6, 3, 5)
e <- 1e-4

g <- function(x) {
  Phi <- pnorm(x) + 1e-10
  res <- dbinom(k, m, Phi, log=TRUE)
  #res <- log(choose(m, k)) + k*log(Phi) + (m-k)*log(1-Phi)
  return(res)
}
dg <- function(x){
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- k*N/Phi - (m-k)*N/(1-Phi)
  return(res)
}
ddg <- function(x){
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- (-k*(N*(x*Phi+N)/Phi^2)) - (m-k)*(N*(N-x*(1-Phi))/(1-Phi)^2)
  return(res)
}
dddg <- function(x){
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  #res <- (-k*(-N*(2*N^2+Phi^2*((f^2)-1)+3*N*x*Phi))/Phi^3) - 
  #       (m-k)*(N*(N*(-2*N+x-x*Phi)+(x^2)*(1-2*Phi+Phi^2) + Phi*(2-Phi)-1)/(1-Phi)^3)
  #res <- (-k*(-2*(N^3)-3*(N^2)*x*Phi-N*(x^2)*(Phi^2)+N*(Phi^2))/(Phi^3)) - 
  #       (m-k)*(N*(2*(N^2)+4*N*x*Phi-N*x*Phi+(x^2)*(Phi^2)-(x^2)+x*N-(Phi^2)+2*Phi-1)/((1-Phi)^3))
  #res <- (-k*(-2*(N^3)-3*(N^2)*x*Phi-N*(x^2)*(Phi^2)+N*(Phi^2))/(Phi^3)) - 
  #       (m-k)*((2*(N^3)-3*(N^2)*x+3*(N^2)*x*Phi+N*(x^2)*(Phi^2)-2*N*(x^2)*Phi-N*(Phi^2)+N*(x^2)+2*N*Phi-N)/((1-Phi)^3))
  #res <- (-k*(-2*(N^3)-3*(N^2)*x*Phi-N*(x^2)*(Phi^2)+N*(Phi^2))/(Phi^3)) - 
  #       (m-k)*(N*(2*(N^2)-3*N*x+3*N*x*Phi+(x^2)*(Phi^2)-2*(x^2)*Phi-(Phi^2)+(x^2)+2*Phi-1)/((1-Phi)^3))
  
  res  <- (-k*(N*(-2*(N^2)-3*N*x*Phi-(Phi^2)*((x^2)-1)))/(Phi^3)) - 
          (m-k)*(N*(N*(2*N-3*x+3*x*Phi)+(x^2)*((Phi^2)-2*Phi+1)-(Phi^2)+2*Phi-1)/((1-Phi)^3))
  return(res)
}
c1 <- checkGrad(g, dg, e, f)
c2 <- checkGrad(dg, ddg, e, f)
c3 <- checkGrad(ddg, dddg, e, f)




#################################
#################################

lp <- function(x) {
  Phi <- pnorm(x) + 1e-10
  return(log(Phi))
}
d1lp  <- function(x) {
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- 1*N / Phi 
  
}
d2lp  <- function(x) {
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- (- N^2 / Phi^2) - (1*x*N)/Phi  
  return(res)
}
d3lp  <- function(x) {
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- (3*x*N^2/Phi^2 + 2*1*N^3/Phi^3 + 1*N*(x^2-1)/Phi)
  return(res)
}
d1 <- checkGrad(lp, d1lp, e, f)
d2 <- checkGrad(d1lp, d2lp, e, f)
d3 <- checkGrad(d2lp, d3lp, e, f)



##############################
###############################
K <- 10
sf2 = 2
lambda=10

se.sf2 <- function(sf2){
  return(sf2^2 * exp(-0.5*K/(lambda^2)))
}
dse.sf2 <- function(sf2){
  return(2 * sf2 * exp(-0.5*K/lambda^2))
}
se.lam <- function(lambda){
  return(sf2^2 * exp(-0.5*K/(lambda^2)))
}
dse.lam <- function(lambda){
  return(sf2^2 * exp(-0.5*K/(lambda^2)) * K/(lambda^3))
}
ds.se <- checkGrad(se.sf2, dse.sf2, e, f)
dl.se <- checkGrad(se.lam, dse.lam, e, f)







a <- .6
b <- .4
c <- .1

f <- c(.4, .5, .7, .2, .9)
t <- c(10, 11, 10, 11, 10)
m <- c(5, 7, 6, 3, 5)
e <- 1e-4

x <- f
LLa <- function(a) {
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  res <- dbinom(m, t, Phi, log=TRUE)
  return(res)
}
L.a <- function(a){
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  N   <- dnorm(g)
  res <- N * x^2 * ( (m-t*Phi)/(Phi*(1-Phi)) )
}
L.aa <- function(a){
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  N   <- dnorm(g)
  res <- x^2 * ( (m * (x^2*N*(-g*Phi - N)/Phi^2) ) - ((t-m)*(x^2*N*(-g+g*Phi+N))/((1-Phi)^2)) )
}


LLb <- function(b) {
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  res <- dbinom(m, t, Phi, log=TRUE)
  return(res)
}
L.b <- function(b){
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  N   <- dnorm(g)
  res <- N * x * ( (m-t*Phi)/(Phi*(1-Phi)) )
}

LLc <- function(c) {
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  res <- dbinom(m, t, Phi, log=TRUE)
  return(res)
}
L.c <- function(c){
  g   <- a*x^2 + b*x + c
  Phi <- pnorm(g) + 1e-10
  N   <- dnorm(g)
  res <- N * ( (m-t*Phi)/(Phi*(1-Phi)) )
}

dl.a  <- checkGrad(LLa, L.a, e, f)
dl.b  <- checkGrad(LLb, L.b, e, f)
dl.c  <- checkGrad(LLc, L.c, e, f)
dll.a  <- checkGrad(L.a, L.aa, e, f)
