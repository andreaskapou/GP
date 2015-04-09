##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
source("checkGrad.R")

f <- c(0, .5, -.7, .2, .9)
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

lp <- function(x, y) {
  Phi <- pnorm(x) + 1e-10
  return(log(Phi))
}
d1lp  <- function(x, y) {
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- 1*N / Phi 
  
}
d2lp  <- function(x, y) {
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- (- N^2 / Phi^2) - (1*x*N)/Phi  
  return(res)
}
d3lp  <- function(x, y) {
  Phi <- pnorm(x) + 1e-10
  N   <- dnorm(x)
  res <- (3*x*N^2/Phi^2 + 2*1*N^3/Phi^3 + 1*N*(x^2-1)/Phi)
  return(res)
}
d1 <- checkGrad(lp, d1lp, e, f)
d2 <- checkGrad(d1lp, d2lp, e, f)
d3 <- checkGrad(d2lp, d3lp, e, f)