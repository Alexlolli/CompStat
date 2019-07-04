#Man and machine collaborations

library(stringi)

#A class to define the insured, and how they behave. Should be rewritten to S4 class.
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



#A method for the insured to generate claims that are eiter fraud or not.
# not finished
GenerateClaims<-function(obj){
  UseMethod("GenerateClaims")
}

GenerateClaims.default <- function(obj) {
  cat("This is a generic that generates a claim simulation for the insured\n")
}

GenerateClaims.Insured <- function(obj){
  obj$Claim<-rbinom(1,1,obj$ClaimProb)
}


#A class for the machine to define how it predicts frauds
Machine<-function(par){
  #Class definition goes here. Preferably S4 class
}


#A class for the operator to define how 
Operator<-function(par){
  #Operator definiton goes here. Preferably S4 class
}