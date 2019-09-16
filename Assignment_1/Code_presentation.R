cat("\014")
rm(list=ls())

library(dplyr)
library(ggplot2)
library(microbenchmark)

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

new_my_Kernel <- function(K, rng = c(-Inf, Inf), BW = "Silverman"){
  if(!is.expression(K)){
    stop("Kernel expression missing",
         call. = FALSE)
  }
  
  obj <- list()
  
  obj$rng = rng

  obj$Kernel <- function(x)
    eval(K)*(x>=rng[1] & x<=rng[2])
  
  Int_K<-try(integrate(obj$Kernel, lower = obj$rng[1], upper = obj$rng[2]), silent = T)
  
  if(class(Int_K)=="try-error"){
    stop(warning(Int_K))
  }else if(!(Int_K$value < 1.01 & Int_K$value > 0.99 & Int_K$abs.error < 0.1)){
    stop(paste("Kernel integrates to", Int_K$value))
  }
  
  obj$DDif <- function(x)
    eval(D(D(K,'x'),'x'))*(x>=rng[1] & x<=rng[2])
  
  obj$sig_2 <- integrate(function(x) x^2* obj$Kernel(x), lower = obj$rng[1], upper = obj$rng[2])
  
  ifelse(obj$sig_2$abs.error < 1e-12,
         obj$sig_2 <- obj$sig_2$value^2, 
         warning("Could not evaluate sigma_K^2"))
  
  obj$K_2 <- integrate(function(x) obj$Kernel(x)^2, lower = obj$rng[1], upper = obj$rng[2])
  
  ifelse(obj$K_2$abs.error < 1e-12,
         obj$K_2 <- obj$K_2$value, 
         warning("Could not evaluate K_2^2"))
  
  obj$silverman <- function(data)
    (4/(3*length(data)))^(1/5)*as.numeric(quantile(data,p=0.75)-quantile(data,p=0.25))/ 1.35
  
  structure(obj , class = c("myKernel"))
}



Ep_kern<-new_my_Kernel(K=expression(3/4*(1-x^2)),rng=c(-1,1), BW = "Silverman")

curve(Ep_kern$Kernel(x),from = -1, to= 1)




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



y<-hist(rnorm(1000),plot=F)

plot(y)

class(y)<-"myKernel"

Bw_select(y)

Bw

Bw_select.myKernel<-function(x)
  print("hej")

plot.histogram



Ep_kern$BW(rnorm(1000))

band<-Ep_kern$BW(x[1:500])

plot(band)


my_Kernel <- function(K,){
  
  #Egenskaber for Kernel
  validate_my_Kernel(new_my_kernel(K,)  )
  
  #BW_selection metode?
}


my_Kernel_Density(kernel, data)


#binning
