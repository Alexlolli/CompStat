cat("\014")
rm(list=ls())
library(dplyr)
library(ggplot2)


all.greater <- function(x, h, ... ){
  i <- FALSE
  for(xx in x){
    if(!(isTRUE(all.equal(xx, h)) | xx < h | is.na(xx) | is.nan(xx))){
      if(i){
        y <- append(y, xx)
      }else{
        i <- TRUE
        y <- xx
      }
    }
  }
  if(!i){
    return(NULL)
  }else{
    return(y)
  }
}

all.greater2 <- function(x, h){
  x[x > h + .Machine$double.eps^0.5 ]
}

all.greater2(seq(0,1,0.1),0.3)

all.greater(c(0.1+0.1+0.1,0.4,5,Inf,-Inf,NA),0.3)

#1.4 ####

infrared <- read.table("C:/Users/Alexander Lollike/Documents/Dropbox/GitHub/CompStat/Data/infrared.txt", header = TRUE)

F12<-infrared$F12

hist(log(F12))
hist(log(F12),breaks = seq(min(log(F12)),max(log(F12)),length.out = 15))
hist(log(F12),breaks = seq(min(log(F12)),max(log(F12)),length.out = sqrt(length(F12))))


# 1.5, 1.6 ####
myBreaks <- function(x, h = 5, type = 1){
  if(!(is.atomic(x) || is.list(x))){
    stop("Input vector is not a vector")
  }
  if(type == 1){
    if(h == 0){
      B <- c(min(x),max(x))
    }else{
      y <- sort(unique(x))
      l <- length(y)
      if( !((l-1)%%h) ){
        B <- c(y[1], sort(y)[(1:(l/h))*h+1]) 
      }else{
        B <- c(y[1], sort(y)[(1:floor((l-1)/h))*h + 1], y[l] )
      }
    }
  }else if(type == 2){
    if(h == 0){
      B <- c(min(x),max(x))
    }else{
      y <- sort(x)
      l <- length(y)
      if( !(l%%h) ){
        B <- unique( c(y[1], y[(1:(l/h))*h ] ))
      }else{
        B <- unique( c(y[1], y[(1:floor((l-h)/h))*h], y[l]))
      }
    }
  }
  B
}

myBreaks(1:11,h=5, type = 2)


#Test
myBreaks(x=c(1,3,2,5,10,11,1,1,3),h=2)
myBreaks(x=1:15,h=2)

hist(log(F12), breaks = myBreaks(log(F12)))
hist(log(F12), breaks = myBreaks(log(F12),type = 2))


myBreaks(x=1:25,h = 5, type = 2)
myBreaks(x=1:26,h = 5, type = 2)

hist(x=1:25, breaks=myBreaks(x=1:25,h = 5, type = 2))
hist(x=1:25, breaks=c(0,5,10,15,20,25))


test <- rnorm(1000)
hist(test, breaks = myBreaks(test, type = 2), freq = T)


#1.7#####
myHist <- function(h = 5, ...){
  hist(log(F12), breaks = myBreaks(log(F12), h = h, type = 2), ...)
}


F12 <- infrared$F12

environment(myHist)

myHist()

rm(F12)

myHist() #fungerer ikke fordi F12 er i global env

F12 <- infrared$F12

myHist()

sum(duplicated(log(F12))) # We should not get feq=5 for all boxes
myHist(h = 5, freq = TRUE) #https://stat.ethz.ch/pipermail/r-help/2008-May/162312.html
#"Basically R is reluctant to let you shoot yourself in the foot unless you are really determined to do so."

myHist(h = 0)

#1.8 ####
environment(myHist)
#global environment

test_envir<-new.env()
environment(myHist) <- test_envir
assign("F12", F12, envir = test_envir)

rm(F12)
myHist()

#1.9 ####
F12 <- infrared$F12

myHist_Factory <- function(x){
  function(h = 5, ...){
    hist(x, breaks = myBreaks(x, h = h, type = 2), ...)
  }
}

myHist_F12 <- myHist_Factory(log(F12))

environment(myHist_F12)

rm(F12)

myHist_F12() #Giver fejl på grund af lazy evaluation! Der er ingen grund til at duplikere F12, det skal påtvinges:
#Enten skal myHist_F12() køres én gang med adgang til F12, eller også skal lazy evaluatio undgås:

F12 <- infrared$F12

myHist_Factory2 <- function(x){
  force(x)
  function(h = 5, ...){
    hist(x, breaks = myBreaks(x, h = h, type = 2), ...)
  }
}

myHist_F12_2 <- myHist_Factory2(log(F12))

environment(myHist_F12_2)

myHist_F12_2()

rm(F12)

myHist_F12_2() #fungerer stadig!


#1.10 ####

tmp <- myHist(10, plot = FALSE)

typeof(tmp) #list
class(tmp) #histogram

plot(tmp, col = "red")
#does Historgram class  have a plot method?
#yes plot.Histogram is a method for objects of class "histogram"

## S3 method for class 'histogram'
#plot(x, freq = equidist, density = NULL, angle = 45,
#     col = NULL, border = par("fg"), lty = NULL,
#     main = paste("Histogram of",
#                  paste(x$xname, collapse = "\n")),
#     sub = NULL, xlab = x$xname, ylab,
#     xlim = range(x$breaks), ylim = NULL,
#     axes = TRUE, labels = FALSE, add = FALSE,
#     ann = TRUE, ...)



tmp <- myHist(10)

typeof(tmp) #list
class(tmp) #histogram
str(tmp)


