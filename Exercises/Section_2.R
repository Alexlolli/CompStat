cat("\014")
rm(list=ls())
library(dplyr)
library(ggplot2)
library(microbenchmark)

infrared <- read.table("C:/Users/Alexander Lollike/Documents/Dropbox/GitHub/CompStat/Data/infrared.txt", header = TRUE)

F12<-infrared$F12

KDens <- function (x, h, K, m = 512) {
  rg <- range(x)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) 
  const <- (h * length(x))
  for (i in seq_along(xx))
    y[i] <- sum(K((x-xx[i])/h)) / const
  list(x = xx, y = y)
}

Gauss_K <- function(x) exp(-x^2/2)/sqrt(2*pi)
Rect_K <- function(x) (x<1 & x>-1)/2
Epan_K <- function(x) 3/4*(1-x^2)*(x<1 & x>-1)

x <- log(F12)

hist(x, prob = T)
lines(KDens(x, h=0.3, K = Gauss_K))
lines(density(x, kernel = "gaussian"))

hist(x, prob = T)
lines(KDens(x, h=0.3, K = Rect_K), lty = 2)
lines(density(x, kernel = "rectangular"))

hist(x, prob = T)
lines(KDens(x, h=0.3, K = Epan_K), lty = 3)
lines(density(x, kernel = "epanechnikov"))

range(KDens(x, h=0.2071, K = Epan_K)$y-density(x, kernel = "epanechnikov")$y)
plot(KDens(x, h=0.2071, K = Epan_K)$y-density(x, kernel = "epanechnikov")$y,type="l")

EpDens <- function (x, h, m = 512) {
  rg <- range(x)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) 
  const <- (h * length(x))
  for (i in seq_along(xx))
    y[i] <- sum(3/4*(1-(x-xx[i])^2/h^2)*(((x-xx[i])/h)<1 & ((x-xx[i])/h)>-1))/const
  list(x = xx, y = y)
}

m<-512
x<-log(F12)
h<-1
rg <- range(x)
xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
y <- numeric(m) 
const <- (h * length(x))
i<-250
z<-(((x-xx[i])/h)<=1 & ((x-xx[i])/h)>=-1)
sum(z)

hist(x, prob = T)
lines(KDens(x, h=0.3, K = Epan_K), lty = 3)
lines(EpDens(x, h=0.3), lty = 3)

tilde_sd <- min( sd(log(F12)) , (quantile(log(F12),p=0.75)-quantile(log(F12),p=0.25)) / 1.35)

f_Gauss_est <- 3/(8*tilde_sd^5*sqrt(pi))

K_2 <- 3/5

sigma_k_4 <- 1/25 

h_n<-(K_2/(sigma_k_4*f_Gauss_est))^(1/5)*length(F12)^(-1/5)

hist(log(F12), prob = T, breaks = 50)
lines(EpDens(log(F12), h = h_n), lty = 3)





x_i<-c(0.2,0.5)
x_j<-c(0.5,0.2)

r<-0.2

pmin(x_j,x_i)-pmax(x_j,x_i)+2*r


x<-log(F12)

silverman <- function(x){
  (4/(3*length(x)))^(1/5)*as.numeric(quantile(x,p=0.75)-quantile(x,p=0.25))/ 1.35
}




f_est_Ep<-function(x,r){
  
  Z <- outer(x,x,function(x_i,x_j){
    z <- pmin(x_j,x_i)-pmax(x_j,x_i)+2*r
    return(z*(z>0))
  })

  S <- sum(Z)*9/4
  S <- S*1/(r^6*length(x)^2)
  S
}


r_1 <- silverman(log(F12))

f_EP_est<-f_est_Ep(log(F12),r_1)


h_n<-(K_2/(sigma_k_4*f_EP_est))^(1/5)*length(F12)^(-1/5)











set.seed(1234)
x <- rnorm(2^13)

conf <- expand.grid(
  fun = c("density","EpDens"),
  k = 2^c(5:13)
)


calls <- paste0(conf[,1], "(x[1:", conf[, 2], "], 0.2)")
expr_list <- lapply(calls, function(x) parse(text = x)[[1]])
kern_benchmarks <- microbenchmark(list = expr_list, times = 40L)

kern_benchmarks

kern_benchmarks

plot(kern_benchmarks)

plot(kern_benchmarks[substr(kern_benchmarks$expr,1,4)=="EpDe",])
plot(kern_benchmarks[substr(kern_benchmarks$expr,1,4)=="dens",], add = T)

