#Man and machine collaborations

library(stringi)

Insured <- function(p) {
  value <- list(name = paste(paste(stri_rand_strings(1, 1, '[A-Z]'),stri_rand_strings(1, 3+floor(5*runif(1)), '[a-z]'),sep=""),
                             paste(stri_rand_strings(1, 1, '[A-Z]'),stri_rand_strings(1, 3+floor(9*runif(1)), '[a-z]'),sep=""), sep=" "),
                predictor = runif(1),
                ClaimProb = p,
                FraudProb = 0.1,
                Claim = NA,
                Fraud = NA)
  # class can be set using class() or attr() function
  class(value) <- "Insured"
  value
}


GenerateClaims<-function(obj){
  UseMethod("GenerateClaims")
}

GenerateClaims.default <- function(obj) {
  cat("This is a generic that generates a claim simulation for the insured\n")
}

GenerateClaims.Insured <- function(obj){
  obj$Claim<-rbinom(1,1,obj$ClaimProb)
}


I<-Insured(0.5)

GenerateClaims(I)

# create list of insured.



rbinom(10,1,0.5)
