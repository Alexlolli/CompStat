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
  

  obj$silverman <- function(data){
    
    n <- length(data)
    
    IQR <- as.numeric(quantile(data,p=0.75)-quantile(data,p=0.25))
    
    sig <- min(sd(data), IQR/1.35)
    
    f_est <- 3/(2*sqrt(pi)*4*sig^5)
    
    (obj$K_2/(obj$sig_2^2*f_est))^(1/5)*n^(-1/5)
  }
  
  obj$BW_plugin <- function(data){
    
    r <- obj$silverman(data)
    
    rng_data <- range(data)
    
    n <- length(data)
    
    s <- numeric(n-1)
    
    int_f<-function(x){
      sapply(x,function(x){
        for(i in seq_along(data[-1])){
          s[i]<-sum(
            obj$DDif((x-data[i])/r)*obj$DDif((x-data[(i+1):n])/r)) #calculate triangle matrix vectorized
        }
        sum(s)*2+sum(obj$DDif((x-data)/r)^2) #Use symmetry of matrix and and diagonal
      })
    }
    
    int_res <- integrate(int_f,lower = rng_data[1]-r, upper = rng_data[2]+r, subdivisions = 100L, rel.tol = 0.05)
    
    f_est <- 1/(n^2*r^6)*int_res$value
    
    (obj$K_2/(obj$sig_2^2*f_est))^(1/5)*n^(-1/5)
  }
  
  obj$LOOCV <- function(data, n, h) {
    
    z <- NULL
    for(i in seq_along(data)) {
      y <- data[-i]
      z[i] <- sum(obj$Kernel((data[i] - y) / h)) / ((n-1) * h)
    }
    sum(z)
  }
  structure(obj , class = c("my_Kernel"))
}

my_Kernel<-function(K, rng = c(-Inf, Inf)){
  new_my_Kernel(K, rng)
}

band_width<-function(data, K, method = "gauss", scale_EP = F){
  
  n <- length(data)

  if(method == "plug_in"){
    
    BW <- K$BW_plugin(data)
    
  }else if(method == "silverman"){
    
    BW <- K$silverman(data)
    
  }else if(method == "LOOCV"){
    
    BW<- optim(par = 5, 
               fn = function(h) - K$LOOCV(data, n, h),
               method = "Brent", lower = 0.0001, upper = 100)$par
    
  }else if(method == "epanechnikov"){
    
    IQR <- as.numeric(quantile(data,p=0.75)-quantile(data,p=0.25))
    
    sig <- min(sd(data),IQR/ 1.35)
    
    f_est <- 3/(2*sqrt(pi)*4*sig^5)
    
    r <- (0.6/(0.04*f_est))^(1/5)*n^(-1/5) #silvermans rule of thumb, epanechnikov pilot
    
    s <- numeric(n-1)
    
    for(i in seq_along(data[-1])){
      z <-  (pmin(data[i], data[(i+1):n])-pmax(data[i], data[(i+1):n])+2*r) #calculate triangle matrix vectorized
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
    IQR <- as.numeric(quantile(data,p=0.75)-quantile(data,p=0.25))
    
    sig <- min(sd(data),IQR/ 1.35)
    
    r <- (4/(3*n))^(1/5)*sig #silvermans rule of thumb, gaussian pilot
    
    s <- numeric(n-1)

    for(i in seq_along(data[-1])){
      #calculate triangle matrix vectorized
      
      z <-  (1 / (2 * sqrt(pi)) * exp(-(data[i] - data[(i+1):n])^2 / (4 * r^2)) *
                                           (3 / 4 * r - 3 / (4 * r) * (data[i] - data[(i+1):n])^2 + 
                                                                      (data[i] - data[(i+1):n])^4 / (16 * r^3)))
      s[i] <- sum(z)
    }
    
    sum(s)
    
    diagonal <- (1 / (2 * sqrt(pi)) * 1 * (3 / 4 * r ))
    
    int_res <- (sum(s)*2+n*diagonal) #Use symmetry of matrix and and diagonal of length 2*r
    
    f_est <- 1/(n^2*r^6)*int_res
    
    if(scale_EP){
      BW <- (K$K_2/(f_est))^(1/5)*n^(-1/5)
    }else{
      BW <- (K$K_2/(K$sig_2^2*f_est))^(1/5)*n^(-1/5)
    }
  }
  BW
}


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
    h <- band_width(x, K, method = bw, scale_EP = scale_EP) 
  }
  
  #First row of toeplitz matrix
  kerneval <- K$Kernel((xx - xx[1])/h)/h
  kerndif <- toeplitz(kerneval)
  y <- colSums(w * kerndif)

  list(x = xx, y = y, h = h)
}



Ep_kern <- my_Kernel(K=expression(3/4*(1-x^2)),rng=c(-1,1))

Ep_kern_sc <- my_Kernel(K=expression(3/4*(1-x^2/5)/sqrt(5)),rng=c(-sqrt(5),sqrt(5)))

G_kern <-my_Kernel(K=expression(1/sqrt(2*pi)*exp(-x^2/2)),rng = c(-Inf, Inf))

BW_GR<-0.9*min(sd(log(F12)),as.numeric(quantile(log(F12),p=0.75)-quantile(log(F12),p=0.25))/1.34)*length(F12)^(-1/5)


dens_G <- my_density(log(F12),G_kern, bw = BW_GR, m=512, scale_EP = F)
dens_GR <- density(log(F12),kernel = "gaussian")
dens_E <- my_density(log(F12),Ep_kern_sc, bw = BW_GR, m=512)
dens_ER <- density(log(F12),kernel = "epanechnikov")

#ingen forskel i bandwidth
dens_GR$bw-BW_GR
dens_ER$bw-BW_GR

#lille forskel i density implementation, muligvis FFT
par(mfrow=c(1,2))
plot(dens_ER$y-dens_E$y)
plot(dens_GR$y-dens_G$y)

range(dens_ER$y-dens_E$y)#epanechnikov
range(dens_GR$y-dens_G$y)#gauss



test_data <- sample(log(F12),512)

BW_GR<-0.9*min(sd(test_data),as.numeric(quantile(test_data,p=0.75)-quantile(test_data,p=0.25))/1.34)*length(test_data)^(-1/5)

dens_G <- my_density(test_data,G_kern, bw = BW_GR, m=512, scale_EP = F)
dens_GR <- density(test_data,kernel = "gaussian")
dens_E <- my_density(test_data,Ep_kern_sc, bw = BW_GR, m=512)
dens_ER <- density(test_data,kernel = "epanechnikov")

#lille forskel i density implementation, muligvis FFT
plot(dens_ER$y-dens_E$y)
plot(dens_GR$y-dens_G$y)

range(dens_ER$y-dens_E$y)#epanechnikov
range(dens_GR$y-dens_G$y)#gauss


dens_G <- my_density_bin(log(F12),G_kern, bw = BW_GR, bins=512)
dens_E <- my_density_bin(log(F12),Ep_kern_sc, bw = BW_GR, bins=512)

plot(dens_ER$y-dens_E$y)
plot(dens_GR$y-dens_G$y)


BW_E <- band_width(log(F12), Ep_kern, method = "epanechnikov")
BW_L <- band_width(log(F12), Ep_kern, method = "LOOCV")
BW_S <- band_width(log(F12), Ep_kern, method = "silverman")
BW_G <- band_width(log(F12), Ep_kern, method = "gauss")

dens_bin_E <- my_density_bin(log(F12),Ep_kern, bw = BW_GR, bins=512, scale_EP = T)
dens_bin_ER <- density(log(F12),kernel = "epanechnikov")

plot(dens_bin_ER$y-dens_bin_E$y)


plot(dens_bin_ER$y-dens_bin_E$y)
hist(log(F12),prob=T, breaks =25)
lines(dens_bin_E$x,dens_bin_E$y)

dens_GR


data <- log(F12)
K<-Ep_kern

w <- 1


BW_E <- band_width(data, Ep_kern, method = "epanechnikov")
BW_L <- band_width(data, Ep_kern, method = "LOOCV")
BW_S <- band_width(data, Ep_kern, method = "silverman")
BW_G <- band_width(data, Ep_kern, method = "gauss")

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





bandwidth_bench <-microbenchmark(band_width(log(F12),Ep_kern,method="silverman"),
                                 band_width(log(F12),Ep_kern,method="epanechnikov"),
                                 band_width(log(F12),Ep_kern,method="gauss"),
                                 band_width(log(F12),Ep_kern,method="LOOCV"))

autoplot(bandwidth_bench) + 
  geom_jitter(position = position_jitter(0.2, 0), 
              aes(color = expr), alpha = 0.4) + 
  aes(fill = I("gray")) + 
  theme(legend.position = "none")


set.seed(1234)
MM_var<-rbinom(Nsim, size = 1,p=0.5) #Bernouilli variables

x <- (1-MM_var)* rnorm(Nsim, mean = par_1[1], sd = par_1[2]) +
  MM_var * rnorm(Nsim, mean = par_2[1], sd = par_2[2])


conf <- expand.grid(
  fun = c("silverman", "epanechnikov", "gauss", "LOOCV"),
  n = 2^(5:12)
)


calls <- paste0('band_width(x[1:', conf[, 2], '], K = Ep_kern, method= "',conf[, 1], '" )', sep = '')
expr_list <- lapply(calls, function(x) parse(text = x)[[1]])
kern_benchmarks <- microbenchmark(list = expr_list, times = 40L)


benchdata<-as.data.frame(cbind(conf,time=summary(kern_benchmarks)$median))


ggplot(benchdata, aes(x = n, y = time, color = fun)) + 
  geom_abline(intercept = 7.5, slope = 1, color = "gray", linetype = 2) +
  stat_summary(fun.y = "median", geom = "line") + 
  stat_summary(fun.y = "median", geom = "point") + 
  facet_wrap(~ "Bandwidth calculation time") + 
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous("Time (s)", trans = "log2", 
                     breaks = c(1e3, 1e4, 1e5, 1e6, 1e7, 1e8), 
                     labels = c("0.001", "0.01", "0.1", "1", "10", "100")) +
  scale_color_discrete("Method:") + 
  theme(legend.position="top")


















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
  obj$BW <- function(data){
    
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
