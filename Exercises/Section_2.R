cat("\014")
rm(list=ls())
library(dplyr)
library(ggplot2)
library(microbenchmark)

#load data

infrared <- read.table("C:/Users/Alexander Lollike/Documents/Dropbox/GitHub/CompStat/Data/infrared.txt", header = TRUE)

F12<-infrared$F12

#general implementation

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

#Epanechnikov kernel

EpDens <- function (x, h, m = 512) {
  rg <- range(x)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) 
  const <- (h * length(x))
  for (i in seq_along(xx))
    y[i] <- sum(3/4*(1-(x-xx[i])^2/h^2)*(((x-xx[i])/h)<1 & ((x-xx[i])/h)>-1))/const
  list(x = xx, y = y)
}


hist(x, prob = T)
lines(KDens(x, h=0.3, K = Epan_K), lty = 3)
lines(EpDens(x, h=0.3), lty = 3)


tilde_sd <- min( sd(log(F12)) , (quantile(log(F12),p=0.75)-quantile(log(F12),p=0.25)) / 1.35)

#Estimation of f_2, assuming true gaussian distribution
f_Gauss_est <- 3/(8*tilde_sd^5*sqrt(pi))

#derived analytically
K_2 <- 3/5
sigma_k_4 <- 1/25 

# Simple bandwidth estimation using weird combination of kernels?
h_n<-(K_2/(sigma_k_4*f_Gauss_est))^(1/5)*length(F12)^(-1/5)

hist(log(F12), prob = T, breaks = 50)
lines(EpDens(log(F12), h = h_n), lty = 3)




#Silvermans rule of thumb. Bandwidth
silverman <- function(x){
  (4/(3*length(x)))^(1/5)*as.numeric(quantile(x,p=0.75)-quantile(x,p=0.25))/ 1.35
}


# Test of silverman

set.seed(1234)
Nsim <- 2^20
MM_var<-rbinom(Nsim, size = 1,p=0.5)
x_1 <- (1-MM_var)*rnorm(Nsim) + MM_var*rnorm(Nsim, mean = 6, sd =2)

hist(x_1)

set.seed(1234)
Nsim <- 2^20
par_1 <- c( 0, 1)
par_2 <- c( 6, 2)


MM_var<-rbinom(Nsim, size = 1,p=0.5) #Bernouilli variables

x <- (1-MM_var)* rnorm(Nsim, mean = par_1[1], sd = par_1[2]) +
  MM_var * rnorm(Nsim, mean = par_2[1], sd = par_2[2])
hist(x, prob = TRUE, breaks = 50, main = paste("Histogram of x. par_1=c(",par_1[1],",",par_1[2],") and par_2=c(",par_2[1],",",par_2[2],")", sep=""))
rg_x <- range(x)
x_vals <- seq(rg_x[1], rg_x[2], length.out = 512)

lines(x=x_vals,y=0.5*dnorm(x_vals, mean = par_1[1], sd = par_1[2])+0.5*dnorm(x_vals, mean = par_2[1], sd = par_2[2]))


x <- lapply(10:20,function(i) x_1[1:2^i])

ans <- lapply(x,function(x) KDens(x, h=silverman(x), K = Gauss_K))

plot(ans[[1]]$y-dnorm(ans[[1]]$x,sd=1,mean=0),type="l")
for(i in c(5,10))
 lines(ans[[i]]$y-dnorm(ans[[i]]$x,sd=1,mean=0),type="l", lty = i)


hist(x_1, prob = T)
lines(ans[[1]])


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



# Niels' nice speed plot


```{r kern-bench-grid, dependson=c("kernDens", "kernDens-vec", "kernDens-apply", "kernDens-outer"), echo=FALSE}
conf <- expand.grid(
  fun = c("kernDens", "kernDens_vec", "kernDens_apply", "kernDens_outer"),
  n = 2^(5:11), 
  m = 2^c(5, 7, 9, 11)
)
set.seed(1234)
x <- rnorm(2^11)
calls <- paste0(conf[, 1], "(x[1:", conf[, 2], "], h = 0.2, m = ", conf[, 3], ")")
expr_list <- lapply(calls, function(x) parse(text = x)[[1]])
kern_benchmarks <- microbenchmark(list = expr_list, times = 40L)
```

```{r kern-bench-fig, dependson="kern-bench-grid", fig.cap="Median run times for the four different implementations of kernel density estimation. The dashed gray line is a reference line with slope 1.", echo=FALSE, warning=FALSE, message=FALSE, out.width="100%"}
class(kern_benchmarks) <- "data.frame"
kern_benchmarks <- bind_cols(conf, expr = calls) %>% 
  left_join(kern_benchmarks, .)
kern_benchmarks$m <- factor(
  kern_benchmarks$m, 
  levels = c(32, 128, 512, 2048),
  labels = c("m = 32", "m = 128", "m = 512", "m = 2048")
)
ggplot(kern_benchmarks, aes(x = n, y = time, color = fun)) + 
  geom_abline(intercept = 15, slope = 1, color = "gray", linetype = 2) +
  stat_summary(fun.y = "median", geom = "line") + 
  stat_summary(fun.y = "median", geom = "point") + 
  facet_wrap(~ m) + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous("Time (ms)", trans = "log2", 
                     breaks = c(1e5, 1e6, 1e7, 1e8), 
                     labels = c("0.1", "1", "10", "100")) +
  scale_color_discrete("Function:") + 
  theme(legend.position="top")
```

