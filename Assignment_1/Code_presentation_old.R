cat("\014")
rm(list=ls())

library(dplyr)
library(ggplot2)
library(microbenchmark)

infrared <- read.table("C:/Users/Alexander Lollike/Documents/Dropbox/GitHub/CompStat/Data/infrared.txt", header = TRUE)

F12<-infrared$F12

set.seed(1234)
Nsim <- 2^20
par_1 <- c( 0, 1)
par_2 <- c( 6, 2)


MM_var<-rbinom(Nsim, size = 1,p=0.5) #Bernouilli variables

x <- (1-MM_var)* rnorm(Nsim, mean = par_1[1], sd = par_1[2]) +
  MM_var * rnorm(Nsim, mean = par_2[1], sd = par_2[2])


rg_x <- range(x)
x_vals <- seq(rg_x[1], rg_x[2], length.out = 512)


hist(x, prob = TRUE, breaks = 60, main = paste("Histogram of x. par_1=c(",par_1[1],",",par_1[2],") and par_2=c(",par_2[1],",",par_2[2],")", sep=""))
lines(x=x_vals,y=0.5*dnorm(x_vals, mean = par_1[1], sd = par_1[2])+0.5*dnorm(x_vals, mean = par_2[1], sd = par_2[2]))


EpDens <- function (x, h, m = 512) {
  rg <- range(x)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) 
  const <- (h * length(x))
  for (i in seq_along(xx))
    y[i] <- sum((1-(x-xx[i])^2/h^2)*(((x-xx[i])/h)<1 & ((x-xx[i])/h)>-1))/const*3/4
  list(x = xx, y = y)
}


#object orient

new_my_Kernel <- function(K, rng){
  if(!is.expression(K)){
    stop("Kernel expression missing",
         call. = FALSE)
  }
  
  obj <- list()
  
  obj$rng = rng

  obj$Kernel <- function(x) #Not idiot proof. Expression in x is needed
    eval(K)*(x>=rng[1] & x<=rng[2])
    
  Int_K<-try(integrate(obj$Kernel, lower = obj$rng[1], upper = obj$rng[2]), silent = T) #Does it integrate to 1?
  
  if(class(Int_K)=="try-error"){
    stop(warning(Int_K))
  }else if(!(Int_K$value < 1.01 & Int_K$value > 0.99 & Int_K$abs.error < 0.1)){
    stop(paste("Kernel integrates to", Int_K$value))
  }
  
  obj$DDif <- function(x) #Double differentiated
    eval(D(D(K,'x'),'x'))*(x>=rng[1] & x<=rng[2])
  
  obj$sig_2 <- integrate(function(x) x^2* obj$Kernel(x), lower = obj$rng[1], upper = obj$rng[2]) #calculate sigma_K^2
  
  obj$K_2 <- integrate(function(x) obj$Kernel(x)^2, lower = obj$rng[1], upper = obj$rng[2]) #Calculate K_2^2
  
  ifelse(obj$sig_2$abs.error < 1e-12,
         obj$sig_2 <- obj$sig_2$value, 
         warning("Could not evaluate sigma_K^2"))
  
  ifelse(obj$K_2$abs.error < 1e-12,
         obj$K_2 <- obj$K_2$value, 
         warning("Could not evaluate K_2^2"))
  

  obj$silverman <- function(data, weights){
    
    
    IQR <- as.numeric(quantile(rep(data,weights),p=0.75)-quantile(rep(data,weights),p=0.25))
    
    sig <- min(sd(rep(data,weights)), IQR/1.35)
    
    f_est <- 3/(2*sqrt(pi)*4*sig^5)
    
    (obj$K_2/(obj$sig_2^2*f_est))^(1/5)*sum(weights)^(-1/5)
  }
  
  obj$BW_plugin <- function(data, weights){
    #Silvermans rule of thumb, adjusted for weights, but not using the actual Kernel! Normal approximation
    r <- obj$silverman(data, weights)
    
    if(length(weights)==1)
      weights<-rep(weights,length(data))
    
    rng_data <- range(data)
    
    n <- sum(weights)
    m <- length(data)
    
    s <- NULL
    
    int_f<-function(x){
      sapply(x,function(x){
        for(i in seq_along(data[-1])){
          s[i]<-sum(
            obj$DDif((x-data[i])/r)*weights[i]*obj$DDif((x-data[(i+1):m])/r)*weights[(i+1):m]) #calculate triangle matrix vectorized
        }
        sum(s)*2+sum(obj$DDif((x-data)/r)^2*weights^2) #Use symmetry of matrix and and diagonal
      })
    }
    
    int_res <- integrate(int_f,lower = rng_data[1]-r, upper = rng_data[2]+r, subdivisions = 100L, rel.tol = 0.05)
    
    f_est <- 1/(n^2*r^6)*int_res$value
    
    (obj$K_2/(obj$sig_2^2*f_est))^(1/5)*n^(-1/5)
  }
  
  obj$LOOCV <- function(data, weights, h) {
    
    if(length(weights)==1)
      weights<-rep(weights,length(data))
    
    z <- NULL
    for(i in seq_along(data)) {
      y <- data[-i]
      w <- weights[-i]
      z[i] <- weights[i]*sum(obj$Kernel((data[i] - y) / h)*w) / (sum(w) * h)
    }
    sum(z)
  }
  structure(obj , class = c("my_Kernel"))
}

my_Kernel<-function(K, rng = c(-Inf, Inf)){
  new_my_Kernel(K, rng)
}

band_width<-function(data, K, weights = 1, method = "gauss", scale_EP = F){
  if(length(weights)==1)
    weights<-rep(weights,length(data))
  
  n <- sum(weights)

  if(method == "plug_in"){
    
    BW <- K$BW_plugin(data, weights)
    
  }else if(method == "silverman"){
    
    BW <- K$silverman(data, weights)
    
  }else if(method == "LOOCV"){
    
    BW<- optim(par = 5, 
               fn = function(h) - K$LOOCV(data, weights, h),
               method = "Brent", lower = 0.0001, upper = 100)$par
    
  }else if(method == "epanechnikov"){
    m <- length(data)
    
    IQR <- as.numeric(quantile(rep(data,weights),p=0.75)-quantile(rep(data,weights),p=0.25))
    
    sig <- min(sd(rep(data,weights)),IQR/ 1.35)
    
    f_est <- 3/(2*sqrt(pi)*4*sig^5)
    
    r <- (0.6/(0.04*f_est))^(1/5)*n^(-1/5) #silvermans rule of thumb, epanechnikov pilot
    
    s <- numeric(m-1)
    
    for(i in seq_along(data[-1])){
      z <- weights[i] * weights[(i+1):m] * (pmin(data[i], data[(i+1):m])-pmax(data[i], data[(i+1):m])+2*r) #calculate triangle matrix vectorized
      s[i] <- sum(z*(z>0))
    }
    
    int_res <- 9/4*(sum(s)*2+n*2*r) #Use symmetry of matrix and and diagonal of length 2*r
    
    f_est <- 1/(n^2*r^6)*int_res
    if(scale_EP){
      BW <- (K$K_2/(K$sig_2^2*f_est))^(1/5)*n^(-1/5)/sqrt(5)
    }else{
      BW <- (K$K_2/(K$sig_2^2*f_est))^(1/5)*n^(-1/5)
    }
    
    
  }else if(method == "gauss"){
    m <- length(data)
    
    IQR <- as.numeric(quantile(rep(data,weights),p=0.75)-quantile(rep(data,weights),p=0.25))
    
    sig <- min(sd(rep(data,weights)),IQR/ 1.35)
    
    r <- (4/(3*n))^(1/5)*sig #silvermans rule of thumb, gaussian pilot
    
    s <- numeric(m-1)

    for(i in seq_along(data[-1])){
      #calculate triangle matrix vectorized
      
      z <-  (1 / (2 * sqrt(pi)) * exp(-(data[i] - data[(i+1):m])^2 / (4 * r^2)) *
                                           (3 / 4 * r - 3 / (4 * r) * (data[i] - data[(i+1):m])^2 + 
                                                                      (data[i] - data[(i+1):m])^4 / (16 * r^3)))
      z <- weights[i] * z
      s[i] <- sum(z)
    }
    
    sum(s)
    
    diagonal <- (1 / (2 * sqrt(pi)) * 1 * (3 / 4 * r ))
    
    int_res <- (sum(s)*2+n*diagonal) #Use symmetry of matrix and and diagonal of length 2*r
    
    f_est <- 1/(n^2*r^6)*int_res
    
    if(scale_EP){
      BW <- (K$K_2/(K$sig_2^2*f_est))^(1/5)*n^(-1/5)/sqrt(5)
    }else{
      BW <- (K$K_2/(K$sig_2^2*f_est))^(1/5)*n^(-1/5)
    }
  }
  BW
}



Ep_kern<-my_Kernel(K=expression(3/4*(1-x^2)),rng=c(-1,1))




kernbin <- function(x, lo, hi, m) {
  w <- numeric(m)
  delta <- (hi - lo) / (m - 1)
  for(i in seq_along(x)) {
    ii <- floor((x[i] - lo) / delta + 0.5) + 1 #we ought to adjust for a-symmetric kernels, but it is not done here
    w[ii] <- w[ii] + 1
  }
  w/sum(w)
}


my_density<-function(x, K, bw = "gauss", m = 512, scale_EP = F){
  rg <- range(x)

  if(is.numeric(bw)){
    h <- bw
  }else{
    h <- band_width(x, K, method = bw, scale_EP = scale_EP) 
  }
  
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  y <- numeric(m) 
  const <- (h * length(x))
  for (i in seq_along(xx))
    y[i] <- sum(K$Kernel((x-xx[i])/h))
  list(x = xx, y = y/const, h = h)
}

my_density_bin<-function(x, K, bw = "gauss", bins = 512, scale_EP =F){
  #find bandwidth first
  rg <- range(x)
  
  w <- kernbin(x, rg[1], rg[2], bins)
  xx <- seq(rg[1], rg[2], length.out = bins)
  
  if(is.numeric(bw)){
    h <- bw
  }else{
    h <- band_width(xx, K, weights = w, method = bw, scale_EP = scale_EP) 
  }
  
  #First row of toeplitz matrix
  kerneval <- K$Kernel((xx - xx[1])/h)/h
  kerndif <- toeplitz(kerneval)
  y <- colSums(w * kerndif)

  list(x = xx, y = y, h = h)
}

hist(log(F12),breaks=25, prob =T)

dens_G <- my_density(log(F12),Ep_kern, bw = "gauss", m=512, scale_EP = T)
dens_E <- my_density(log(F12),Ep_kern, bw = "epanechnikov", m=512, scale_EP = T)
dens_GR <- density(log(F12),kernel = "gaussian")
dens_ER <- density(log(F12),kernel = "epanechnikov")

plot(dens_GR$y-dens_G$y)

BW_E <- band_width(log(F12), Ep_kern, weights = 1, method = "epanechnikov")
BW_L <- band_width(log(F12), Ep_kern, weights = 1, method = "LOOCV")
BW_S <- band_width(log(F12), Ep_kern, weights = 1, method = "silverman")
BW_G <- band_width(log(F12), Ep_kern, weights = 1, method = "gauss")

dens_bin_G <- my_density_bin(log(F12),Ep_kern, bw = BW_G, bins=512, scale_EP = T)
dens_bin_E <- my_density_bin(log(F12),Ep_kern, bw = BW_E, bins=512, scale_EP = T)
dens_bin_GR <- density(log(F12),kernel = "gaussian")
dens_bin_ER <- density(log(F12),kernel = "epanechnikov")

plot(dens_bin_GR$y-dens_bin_G$y)


plot(dens_bin_ER$y-dens_bin_E$y)
hist(log(F12),prob=T, breaks =25)
lines(dens_bin_E$x,dens_bin_E$y)

dens_GR


data <- log(F12)
K<-Ep_kern

w <- 1


BW_E <- band_width(data, Ep_kern, weights = w, method = "epanechnikov")
BW_L <- band_width(data, Ep_kern, weights = w, method = "LOOCV")
BW_S <- band_width(data, Ep_kern, weights = w, method = "silverman")
BW_G <- band_width(data, Ep_kern, weights = w, method = "gauss")

BW_E
BW_L
BW_G
BW_S


hist(data,breaks=50, prob = T)
lines(EpDens(data, h = BW_P))
lines(EpDens(data, h = BW_L))
lines(EpDens(data, h = BW_E))


hist(x,breaks=50, prob = T)
lines(EpDens(x, h = BW_P))
lines(EpDens(x, h = BW_L))
lines(EpDens(x, h = BW_E))

























if(BW == "Silverman"){
  obj$BW <- obj$silverman
  
}else if(BW == "Plug_in"){
  obj$BW <- function(data){
    r <- obj$silverman(data)
    
    rng_data <- range(data)
    
    int_f2<-function(x)
      sapply(x,function(x)
        sum(outer(data,data,function(y,z) obj$DDif((x-y)/r)*obj$DDif((x-z)/r))))
    
    int_res <- integrate(int_f2,lower = rng_data[1]-r, upper = rng_data[2]+r, subdivisions = 100L, rel.tol = 0.05)
    
    f_est <- 1/(length(data)^2*r^6)*int_res$value
    
    (obj$K_2/(obj$sig_2^2*f_est))^(1/5)*length(data)^(-1/5)
    
  }
}else if(BW == "Plug_in_2"){
  obj$BW <- function(data, weights = 1){
    
    r_last <- obj$silverman(data)
    
    rng_data <- range(data)
    
    r<-r_last
    
    bw_save<-NULL
    
    bw_save[1]<-r
    
    for(i in 2:16){ # Iterate 15 times
      int_f2<-function(x)
        sapply(x,function(x)
          sum(outer(data,data,function(y,z) obj$DDif((x-y)/r_last)*obj$DDif((x-z)/r_last))))
      
      int_res <- integrate(int_f2,lower = rng_data[1]-r, upper = rng_data[2]+r, subdivisions = 100L, rel.tol =  0.07*0.9^i)
      
      f_est <- 1/(length(data)^2*r^6)*int_res$value
      
      r <- (obj$K_2/(obj$sig_2^2*f_est))^(1/5)*length(data)^(-1/5)
      
      smaller <- 0L
      bigger <- 0L
      
      if((r_last-r)>0){
        smaller <- smaller + 1
        if(smaller == i){# we have not yet found lower bound
          r <- r/2 #make next guess smaller
        }else{
          r <- (r_last-r)/2 #make next guess smaller
        }
      }else{
        bigger  <- bigger + 1
        if(bigger == i){# we have not yet found lower bound
          r <- r*2 #make next guess bigger
        }else{
          r <- (r-r_last)/2 #make next guess bigger
        }
        
      }
      bw_save[i]<-r
    }
    bw_save
  }
}else if(BW == "Cross_validation"){
  
}else{
  warning("No bandwidth selection method given")
}
