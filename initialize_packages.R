initialize_packages <- function(work=NULL) {
  if (is.null(work)) {
    suppressMessages(library(pracma, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(data.table, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(boot, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(parallel, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(dplyr, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(tictoc, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(stringr, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(optimParallel, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(R.utils, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(microbenchmark, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(argparser, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(robustbase, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(foreach, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(doParallel, quietly=TRUE, warn.conflicts=FALSE))
    suppressMessages(library(progress, quietly=TRUE, warn.conflicts=FALSE))

  }
}

# Calling the Function for Command Line
initialize_packages()