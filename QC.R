#' Allelic Specific Expression (ASE) Quality Control Test
#'
#' This test is designed to detect poor quality ASE data. This includes
#' low coverage samples, samples with high amounts of monoallelically expressed
#' genes, samples with too wide distributions or skew.
#'
#' There are a few return parameters. QC will always write a file of the
#' samples to keep.
#'
#' @param ref Dataframe containing reference count data
#' @param alt Dataframe containing alternate count data
#' @param numCores number of cores used for parallelization (default = 1). One sample executed per core.
#' @return Dataframe with fitted standard deviations, and coverage, and the cutoff to classify a sample as good quality.


QC <- function(ref, alt, numCores = 1) {
  # input validation
  if (!is.numeric(numCores)) {
    stop("Numeric parameter contains non numeric type")
  }
  if (!is.data.frame(ref) || !is.data.frame(alt)) {
    stop("Either ref or alt is not of type dataframe")
  }
  if (ncol(ref) != ncol(alt) || nrow(ref) != nrow(alt)) {
    stop("Datasets contain different dimensions")
  }
  # Initializing dataframe
  samples_output_list <- colnames(ref)[-2]

  # Temporary Files
  file_out <- tempfile(pattern = "file.out.txt")
  file_indices <- tempfile(pattern = "file.indices.txt")

  # Remove Temporary Files
  if (file.exists(file_out)) {
    # Delete file if it exists
    file.remove(file_out)
  }

  if (file.exists(file_indices)) {
    # Delete file if it exists
    file.remove(file_indices)
  }

  cluster_for_loop <- makeCluster(numCores)
  registerDoParallel(cluster_for_loop)

  foreach(i = 1:ncol(ref), .packages = "pracma", .combine = rbind) %dopar% {
    source("blnmix.R", local = TRUE)
    source("bln.R", local = TRUE)
    source("helper.R", local = TRUE)
    source("internal.R", local = TRUE)
    source("preprocess.R", local = TRUE)
    source("QC.R", local = TRUE)
    source("RcppExports.R", local = TRUE)

    # Getting Sample
    if (!is.numeric(ref[, i]) || !is.numeric(alt[, i])) {
      writeLines("Non numerical Row: Moving to next Row")
    } else {
      r <- (ref[, i])
      a <- (alt[, i])
      t <- r + a
      lst <- remove_total(r, t)
      params <- Fit_BLN_uniform_mixture(lst$ref, lst$total, numCores = numCores)
      # Saving to vectors
      params <- c(params, lst$num.kept, median(lst$total))

      # Write to file
      write(params, file_out, append = TRUE)
      write(i, file_indices, append = TRUE)
    }
  }
  stopCluster(cluster_for_loop)

  # Reading written output
  df_output <- read.table(file_out, sep = " ", fill=TRUE)

  sample_indices <- read.table(file_indices, sep = " ", fill=TRUE)

  # Setting the file index
  df_output$indices <- sample_indices$V1

  # Sorting in Order of original input
  df_output <- setorder(df_output, indices)


  # Adding in Samples
  df_output$sample <- samples_output_list

  # Changing column names
  colnames(df_output) <- c("mean", "std", "num.kept", "median.coverage", "indices", "sample")

  cols_to_keep <- c("mean", "std", "num.kept", "median.coverage", "sample")

  df_output <- df_output[, cols_to_keep]

  # Remove fit for index row
  df_output <- df_output %>% slice(-1)


  # Delete Temporary Files
  if (file.exists(file_out)) {
    # Delete file if it exists
    file.remove(file_out)
  }

  if (file.exists(file_indices)) {
    # Delete file if it exists
    file.remove(file_indices)
  }

  # Return File
  threshold <- adjboxStats(df_output$std)$stats[5]


  df_output$cutoff <- threshold

  return(df_output)
}
