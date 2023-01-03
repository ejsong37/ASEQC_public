#' Probability Density Function (PDF) for the Binomial Logit-Normal
#' Uniform Mixture Model
#'
#' @param ref_counts Reference read counts.
#' @param total_counts Total read counts.
#' @param mean mean of binomial
#' @param sd standard deviation
#' @param lambda lambda (latent variable)
#' @return probability from the BLNMM
dblnmix <- function(x,size,mean=0,sd=1,lambda=1) {
  BLN_L <- dbln(x, size, mean, sd)
  UNIF_L <- (1/(size + 1))
  PROB <- lambda * BLN_L + (1 - lambda) * UNIF_L
  return(PROB)
}

#' Generates data from a Binomial Logit-Normal Uniform Mixture Model (BLNMM)
#'
#' @param n Number of datapoints
#' @param size Total read counts.
#' @param mean mean of the BLNMM
#' @param sd standard deviation of the BLNMM
#' @param lambda lambda of the BLNMM
#' @return Datapoints generated from the BLNMM
rblnmix <- function(n, size,mean=0,sd=1,lambda=0.5) {
  count <- c()
  dist <- c()
  if(length(size) == 1) {
    for(i in seq(n)) {
      vec <- rblnmix_one(size,mean,sd,lambda)
      dist <- c(dist, vec[1])
      count <- c(count, vec[2])
    }
    return(data.frame(counts=count,distribution=dist))
  } else {
    if(n != length(size)) {
      stop('n and length of size do not match')
    }
    for(i in seq(n)) {
      vec <- rblnmix_one(size[i],mean,sd,lambda)
      dist <- c(dist, vec[1])
      count <- c(count, vec[2])
    }
    return(data.frame(counts=count,total=size,distribution=dist))
  }


}

#' Generates one datapoint from a Binomial Logit-Normal Uniform Mixture Model (BLNMM)
#'
#' @param size Total read counts.
#' @param mean mean of the BLNMM
#' @param sd standard deviation of the BLNMM
#' @param lambda lambda of the BLNMM
#' @return datapoint generated from BLNMM
rblnmix_one <- function(size,mean,sd,lambda) {

  # Coin toss simulation
  sim <- sample(100,1)
  lambda <- lambda*100

  if(sim <= lambda) { # Simulate from BLN
    return(c("BLN",rbln(1,size,mean,sd)))
  } else { # Simulate from Uniform
    return(c("Uniform",sample(0:size+1,1,replace=T)))
  }
}
