
 
bandwidth <- function(data, method = "plug-in", pilot = "gauss") {
  
  if(method == "silverman"){
      norm_K <- 0.6
      sigma_K <- 0.2
      IQR <- as.numeric(quantile(data, p = 0.75) - 
                          quantile(data, p = 0.25))
      sigma <- min(sd(data), IQR / 1.35)
      ddf <- 1 / (2 * sqrt(pi)) * 3 / (4 * sigma^5)
      
      h <- (norm_K / (sigma_K^2 * ddf))^(1 / 5) * 
        length(data)^(- 1 / 5)
    }
    
  if(method == "plug-in") {
    
    n <- length(data)
    s <- numeric(n-1)
    if(pilot == "gauss") {
      #first calculate silverman for gaussian as r
      norm_K <- 1 / (2 * sqrt(pi))
      sigma_K <- 1
      IQR <- as.numeric(quantile(data, p = 0.75) - 
                          quantile(data, p = 0.25))
      sigma <- min(sd(data), IQR / 1.35)
      ddf <- 1 / (2 * sqrt(pi)) * 3 / (4 * sigma^5)
      
      r <- (norm_K / (sigma_K^2 * ddf))^(1 / 5) * 
        length(data)^(- 1 / 5)
      # make plug in estimator of ||f||
      for(i in seq_along(data[-1])){
        y <- data[i] - data[(i + 1):n]
        z <- 1 / (2 * sqrt(pi)) * exp(-y^2 / (4 * r^2))
        zz <- 3 / 4 * r - 3 / (4 * r) * y^2
        zz <- zz + y^4 / (16 * r^3)
        z <- z * zz
        s[i] <- sum(z)
      }
      diagonal <- (1 / (2 * sqrt(pi)) * 1 * (3 / 4 * r))
      f <- 1 / (n^2 * r^6) * (sum(s) * 2 + n * diagonal)
      # optimal bandwidth
      h <- (0.6 / (0.2^2 * f))^(1 / 5) * n^(- 1 / 5)
    }
    
    if(pilot == "ep") {
      #first calculate silverman for Epanechnikov as r
      norm_K <- 0.6
      sigma_K <- 0.2
      IQR <- as.numeric(quantile(data, p=0.75) - 
                          quantile(data, p=0.25))
      sigma <- min(sd(data), IQR/ 1.35)
      ddf <- 1 / (2 * sqrt(pi)) * 3 / (4 * sigma^5)
      
      r <- (norm_K / (sigma_K^2 * ddf))^(1 / 5) * 
        length(data)^(- 1 / 5)
      
      for(i in seq_along(data[-1])){
        z <- pmin(data[i], data[(i + 1):n])
        z <- z - pmax(data[i], data[(i + 1):n])
        z <- z + 2 * r
        s[i] <- sum(z * (z > 0))
      }
      f <- 1 / (n^2 * r^6) * 9 / 4 * 
        (sum(s) * 2 + n * 2 * r)
      
      h <- (norm_K / (sigma_K^2 * f))^(1 / 5) * n^(- 1 / 5)
    }
  }
  
  if(method == "cv") {
    # Epanechnikov kernel
    kern <- function(x)
    (abs(x) <= 1) * (1 - x^2) * 3 / 4
    #make the lCV function 
    LCV <- function(data, h) {
      n <- length(data)
      z <- numeric(n)
      for(i in seq_along(data)) {
        y <- data[-i]
        z[i] <- sum(kern((data[i] - y) / h)) / 
          ((n - 1) * h)
      }
      sum((z))
    }
    #optimize the function for h
    h <- optim(par = 0.2, 
               fn = function(h) - LCV(data, h),
               method = "Brent", lower = 0.000001, 
               upper = 10^3)$par
  }
  
  h
}






