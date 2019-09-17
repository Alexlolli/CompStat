
myKernel <- function(K, rng = c(-Inf, Inf)){
  
  if(!is.expression(K)){
    stop("Kernel expression missing",
         call. = FALSE)
  }
  
  obj <- list()
  
  obj$rng = rng
  
  obj$kernel <- function(x)
    eval(K) * (x >= rng[1] & x <= rng[2])
  
  Int_K <- try(integrate(obj$kernel, lower = obj$rng[1], 
                         upper = obj$rng[2]), silent = T)
  
  if(class(Int_K) == "try-error"){
    stop(warning(Int_K))
  }else if(!(Int_K$value < 1.01 & Int_K$value > 0.99 & 
             Int_K$abs.error < 0.1)){
    stop(paste("Kernel integrates to", Int_K$value))
  }
  
  #obj$DDif <- function(x)
  #  eval(D(D(K, 'x'), 'x'))*(x >= rng[1] & x <= rng[2])
  
  obj$sig_2 <- integrate(function(x) x^2 * obj$kernel(x), 
                         lower = obj$rng[1], 
                         upper = obj$rng[2])
  
  ifelse(obj$sig_2$abs.error < 1e-5,
         obj$sig_2 <- obj$sig_2$value, 
         warning("Could not evaluate sigma_K^2"))
  
  obj$K_2 <- integrate(function(x) obj$kernel(x)^2, 
                       lower = obj$rng[1], 
                       upper = obj$rng[2])
  
  ifelse(obj$K_2$abs.error < 1e-5,
         obj$K_2 <- obj$K_2$value, 
         warning("Could not evaluate K_2^2"))
  
  obj$silverman <- function(data, norm_K, sigma_K) {
    sigma <- min(sd(data), 
                 as.numeric(quantile(data, p=0.75) - 
                                        quantile(data, p=0.25))/ 1.35)
    ddf <- 1 / (2 * sqrt(pi)) * 3 / (4 * sigma^5)
    (norm_K / (sigma_K^2 * ddf))^(1 / 5) * 
      length(data)^(- 1 / 5)
  }
  
  obj$LOOCV <- function(data, h) {
    n <- length(data)
    z <- numeric(n)
    for(i in seq_along(data)) {
      y <- data[-i]
      z[i] <- sum(obj$kernel((data[i] - y) / h)) / 
        ((n - 1) * h)
    }
    sum((z))
  }
  
  structure(obj , class = c("myKernel"))
}

bandwidth.myKernel <- function(obj, data, 
                               method = "plug-in",
                               pilot = "gauss") {
  
  method <- match.arg(method, c("silverman", "plug-in", 
                                "cv"))
  # Fejlmeddelelse hvis metode ikke findes
  
  # Tjek at obj er af class myKernel
  
  if(method == "silverman")
    h <- obj$silverman(data, obj$K_2, obj$sig_2)
  
  if(method == "plug-in") {
    
    n <- length(data)
    s <- numeric(n-1)
    if(pilot == "gauss") {
      r <- obj$silverman(data, 1 / (2 * sqrt(pi)), 1)
      for(i in seq_along(data[-1])){
                z <- (1 / (2 * sqrt(pi)) * 
                 exp(-(data[i] - data[(i + 1):n])^2 / 
                       (4 * r^2)) * (3 / 4 * r - 3 / (4 * r) * 
                                       (data[i] - data[(i + 1):n])^2 +
                    (data[i] - data[(i + 1):n])^4 / (16 * r^3)))
        s[i] <- sum(z)
      }
      diagonal <- (1 / (2 * sqrt(pi)) * 1 * (3 / 4 * r ))
      f <- 1 / (n^2 * r^6) * (sum(s) * 2 + n * diagonal)
      h <- (obj$K_2 / (obj$sig_2^2 * f))^(1 / 5) * n^(- 1 / 5)
    }
    
    if(pilot == "ep") {
      r <- obj$silverman(data, 0.6, 0.2)
      for(i in seq_along(data[-1])){
        z <- (pmin(data[i], data[(i + 1):n]) - 
                pmax(data[i], data[(i + 1):n]) + 2 * r)
        s[i] <- sum(z * (z > 0))
      }
      f <- 1 / (n^2 * r^6) * 9 / 4 * 
        (sum(s) * 2 + n * 2 * r)
        
      h <- (obj$K_2 / (obj$sig_2^2 * f))^(1 / 5) * n^(- 1 / 5)
    }
  }
  
  if(method == "cv") {
    h <- optim(par = 0.2, 
               fn = function(h) - obj$LOOCV(data, h),
               method = "Brent", lower = 0.000001, 
               upper = 10^3)$par
  }
  
  h
}

kernbin <- function(x, lo, hi, m) {
  w <- numeric(m)
  delta <- (hi - lo) / (m - 1)
  for(i in seq_along(x)) {
    ii <- floor((x[i] - lo) / delta + 0.5) + 1
    w[ii] <- w[ii] + 1
  }
  w / sum(w)
}


density.myKernel <- function(obj, data, m = 512,
                             bw = "plug-in",
                             pilot = "gauss",
                             binning = FALSE) {
  
  if(is.numeric(bw)) h <- bw
  else h <- bandwidth.myKernel(obj, data, method = bw, pilot)
  
  rg <- range(data) + c(-3 * h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = m)
  
  if(binning == FALSE){
    y <- numeric(m) 
    for (i in seq_along(xx)) {
      y[i] <- sum(obj$kernel((xx[i] - data) / h))
    }
    y <- y / (length(data) * h)
  }
  if (!binning == FALSE){
    weights <- kernbin(data, rg[1], rg[2], m)
    kerneval <- obj$kernel((xx - xx[1]) / h) / h
    kerndif <- toeplitz(kerneval)
    y <- colSums(weights * kerndif)
  }
  list(x = xx, y = y, h = h)
}


