#' Fits Binomial Logit-Normal Uniform Mixture Model (BLNMM)
#'
#' @param ref_counts Reference read counts.
#' @param total_counts Total read counts.
#' @param start Optional Starting point for optimization (Mean, std, lambda)
#' @param factr relative tolerance
#' @param numCores Number of cores for parallelization
#' @return Estimated parameters generated from the fit



Fit_BLN_uniform_mixture <- function(ref, total, start,min_coverage=5,max_coverage=5000,factr=1e7,numCores = 1) {
  if(missing(start)) {
    start.Std <- 1 # starting point for sigma prediction
    start.Mean <- 0 # Starting point for mu prediction
    start.Lambda <- 0.75 # Starting point for lambda
  } else {
    start.Mean <- start[1]
    start.Std <- start[2]
    start.Lambda <- start[3]
  }
  
  print(length(ref))
  # totalcores <- detectCores()
  # stopifnot("Number of cores must not exceed total cores!" = numCores <= totalcores)
  # cl <- makeCluster(numCores)
  # print('Setting up Cluster')
  # printf("Using %d clusters\n",numCores)
  # setDefaultCluster(cl=cl)
  # clusterExport(cl,list('t', 'BLN_uniform_mixture_likelihood',
                        # "dbln","max.length","logit","eps","logistic",
                        # "pln","erf","fxpdf",'remove_total','beep'))

  # MLE_optim <- optimParallel(par <- c(start.Mean,start.Std), fn = BLN_uniform_mixture_likelihood2, ref = ref, total = total,lower = c(-0.25,0), upper = c(0.25,Inf), method = "L-BFGS-B", control=list(factr = factr), parallel=c(cl=cl))
  MLE_optim <- optim(par <- c(start.Mean,start.Std), fn = BLN_uniform_mixture_likelihood2, ref = ref, total = total,lower = c(-0.25,0), upper = c(0.25,Inf), method = "L-BFGS-B")  # Optimizer function
  
print(MLE_optim)
 
#  beep()

# setDefaultCluster(cl=NULL)
# stopCluster(cl)

  return(MLE_optim$par)
}

#' Likelihood function for the Binomial Logit-Normal Uniform Mixture Model (BLNMM)
#'
#' @param ref Reference read counts.
#' @param total Total read counts.
#' @param params Parameters for the mixture model (mean, std, lambda)
#' @return Likelihood calculated for the given parameters.
BLN_uniform_mixture_likelihood <- function(params,ref, total) {
  # parameters
  mean <- params[1]
  sd <- params[2]
  lambda <- params[3]


  # Calculating Likelihood
  BLN_L <- dbln(ref, total, mean, sd)
  UNIF_L <- (1/(total + 1))
  LIKE <- lambda * BLN_L + (1 - lambda) * UNIF_L
  NLL <- -sum(log(LIKE))
  return(NLL)
}

#' Likelihood function for the Binomial Logit-Normal Uniform Mixture Model (BLNMM)
#'
#' @param ref Reference read counts.
#' @param total Total read counts.
#' @param params Parameters for the mixture model (mean, std, lambda)
#' @return Likelihood calculated for the given parameters.
BLN_uniform_mixture_likelihood2 <- function(params,ref, total) {
  # parameters
  mean <- params[1]
  sd <- params[2]


  # Calculating Likelihood
  BLN_L <- dbln(ref, total, mean, sd)
  UNIF_L <- (1/(total + 1))
  LIKE <- .999 * BLN_L + (1 - .999) * UNIF_L
  NLL <- -sum(log(LIKE))
  return(NLL)
}

#' Likelihood function for the Binomial Logit-Normal Distribution
#'
#' @param ref Reference read counts.
#' @param total Total read counts.
#' @param params Parameters for the mixture model (mean, std)
#' @return Likelihood calculated for the given parameters.
BLN_likelihood <- function(params,ref, total) {
  # parameters
  mean <- params[1]
  sd <- params[2]

  # Calculating Likelihood
  print(mean)
  print("-")
  print(sd)
  print('--')
  LIKE <- dbln(ref, total, mean, sd)

  NLL <- -sum(log(LIKE))
  return(NLL)
}

#' Likelihood function for the Binomial Logit-Normal Distribution
#'
#' @param ref Reference read counts.
#' @param total Total read counts.
#' @param params Parameters for the mixture model (mean, std)
#' @return Likelihood calculated for the given parameters.
BLN_Mixture_likelihood <- function(params,ref, total) {
  # parameters
  mean <- params[1]
  sd <- params[2]
  lambda <- params[3]
  print(mean)
  print("-")
  print(sd)
  print('--')
  # Calculating Likelihood
  BLN_L <- dbln(ref, total, mean, sd)
  BLN2_L <- dbln(ref, total, 8, sd) + dbln(ref, total, -8, sd)
  LIKE <- lambda * BLN_L + ((1 - lambda)/2) * BLN2_L
  NLL <- -sum(log(LIKE))
  return(NLL)
}



#' Interquartile Quartile Range (IQR) test for outliers used in ASE_QC.
#' Any value outside the bounds is classified as an outlier.
#' Bounds are calculated as follows (Q1 - 1.5IQR, Q3 + 1.5IQR)
#'
#' @param df Dataframe of fits from ASE_QC
#' @param file Optional parameter for name of the output file
#' @param bounds Whether to write to a file the outlier bounds for each parameter
#' @return dataframe of fits within the outlier bounds
IQR_test <- function(df,file,bounds=F) {

  if(missing(file)) {
    warning("No file name given. Using good_samples.tsv instead")
    file <- "good_samples.tsv"
  }
  # mean test
  mean_bounds <- c(quantile(df$mean, 0.25) - 1.5*IQR(df$mean), quantile(df$mean,0.75) + 1.5*IQR(df$mean))
  # std test
  std_bounds <- c(quantile(df$std, 0.25) - 1.5*IQR(df$std), quantile(df$std,0.75) + 1.5*IQR(df$std))
  # Number of genes kept test
  num_bounds <- c(quantile(df$num.kept, 0.25) - 1.5*IQR(df$num.kept), quantile(df$num.kept,0.75) + 1.5*IQR(df$num.kept))
  # coverage test
  coverage_bounds <- c(quantile(df$median.coverage, 0.25) - 1.5*IQR(df$median.coverage), quantile(df$median.coverage,0.75) + 1.5*IQR(df$median.coverage))
  # lambda test
  lambda_bounds <- c(quantile(df$lambda, 0.25) - 1.5*IQR(df$lambda), quantile(df$lambda,0.75) + 1.5*IQR(df$lambda))

  df <- df[which((df$mean > mean_bounds[1] & df$mean < mean_bounds[2]) &
                   (df$std > std_bounds[1] & df$std < std_bounds[2]) &
             (df$num.kept > num_bounds[1] & df$num.kept < num_bounds[2]) &
             (df$median.coverage > coverage_bounds[1]
              & df$median.coverage < coverage_bounds[2]) &
             (df$lambda > lambda_bounds[1] & df$lambda < lambda_bounds[2])),]

  df2 <- data.frame(mean=mean_bounds,std=std_bounds,num=num_bounds,coverage=coverage_bounds,lambda=lambda_bounds)
  if(bounds==F) {
    write.table(df2,file='bounds.tsv',sep="\t")
  }
  write.table(df, file = file, row.names=FALSE, sep="\t")
  return(df)
}

