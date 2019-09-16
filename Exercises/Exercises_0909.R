cat("\014")
rm(list=ls())
library(dplyr)
library(ggplot2)

infrared <- read.table("C:/Users/Alexander Lollike/Documents/Dropbox/GitHub/CompStat/Data/infrared.txt", header = TRUE)

F12<-infrared$F12

myBreaks <- function(x, h = 5) {
  x <- sort(x)
  breaks <- xb <- x[1]
  k <- 1
  for(i in seq_along(x)[-1]) {
    if (k < h) {
      k <- k + 1
    } else {
      if (xb < x[i - 1] && x[i - 1] < x[i]) {
        xb <- x[i - 1]
        breaks <- c(breaks, xb)
        k <- 1
      }
    }
  }
  last <- length(breaks)
  if(k == min(h, length(x) - 1)) last <- last + 1
  breaks[last] <- x[length(x)]
  breaks
}

myHist_old <- function(h = 5, ...){
  hist(log(F12), breaks = myBreaks(log(F12), h = h), ...)
}


myHist <- function(h = 5, ...){
  x <- hist(log(F12), breaks = myBreaks(log(F12), h = h), ...)
  x <- structure(x , class = c("myHistogram"))
  invisible(x)
}

myHist()

print.myHistogram <- function(x){
  cat("Number of breaks =", length(x$breaks))
}

print(myHist())

plot(myHist())


myHist <- function(h = 5, ...){
  x <- hist(log(F12), breaks = myBreaks(log(F12), h = h), ...)
  x$bins <- length(x$mids)
  invisible(structure(x , class = c("myHistogram", "histogram")))
}

plot(myHist())

summary.myHistogram <- function(x){
  data.frame(mids = x$mids, counts = x$counts)  
}

summary(myHist())

plot.myHistogram <- function(x, ...){
  d <- data.frame(x1=x$breaks[-length(x$breaks)], x2=x$breaks[-1], y1=rep(0,x$bins), y2=x$density, t=rep('a', x$bins), r=1:x$bins)
  p <- geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), ... )
  ggplot() + 
    labs(title = paste("Histogram of log(F12)"),
         y = "Density",
         x = "x") +
    theme(legend.position = "none") +
    p
}


plot(myHist(40), alpha = 0.5)
