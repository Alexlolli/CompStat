##### Test implementation #####

myKernel_ep <- myKernel(K = expression(3 / 4 * (1 - x^2 / 5) / sqrt(5)),
                        rng = c(-sqrt(5), sqrt(5)))

myKernel_gauss <- myKernel(K = expression(1 / (sqrt(2 * pi)) * exp(-x^2 / 2)),
                           rng = c(-Inf, Inf))

infrared <- read.table("infrared.txt", header = TRUE)
F12 <- infrared$F12

bw_dens <- density(log(F12), kernel = "epanechnikov")$bw
bandwidth.myKernel(myKernel_ep, log(F12), method = "silverman") 

bw_dens <- density(log(F12), kernel = "gaussian")$bw

bandwidth.myKernel(myKernel_gauss, log(F12), method = "silverman") 


range(density(log(F12), kernel = "epanechnikov")$y - 
        density.myKernel(myKernel_ep, log(F12), bw = bw_dens)$y)

par(mfrow = c(1, 2))

plot(density(log(F12), kernel = "epanechnikov")$y, type = "l", col = "red",
     xlab = "x", ylab = "Density", lwd = 4)
lines(density.myKernel(myKernel_ep, log(F12), bw = bw_dens)$y, 
      col = "blue", lwd = 2)

plot(density(log(F12), kernel = "epanechnikov")$y - 
       density.myKernel(myKernel_ep, log(F12), bw = bw_dens)$y, 
     type = "l", lwd = 1, 
     ylab = "Difference", xlab = "x", main = "")

par(mfrow = c(1, 1))
