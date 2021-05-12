cuslogit <- function(x){return(log(x/(1-x)))} #transforms between 0 and 1, but can still return NAN's
invcuslogit <- function(x) {return( 1/(1+exp(-x)))} #transforms between 2 and 0  
