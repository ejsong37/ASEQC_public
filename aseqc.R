# This is a wrapper script to take in an input dataset and run ASEQC

# Install Packages 
source("initialize_packages.R")

# Load packages
source("blnmix.R")
source("bln.R")
source("helper.R")
source("internal.R")
source("preprocess.R")
source("QC.R")
source("RcppExports.R")

# Initialize Parser
parser <- arg_parser("Input Data for ASEQC")
num_cores <- detectCores()

parser <- add_argument(parser, "--refcounts", help="Data (TSV) of Reference Counts for Genes")
parser <- add_argument(parser, "--altcounts", help="Data (TSV) of Alternate Counts for Genes")
parser <- add_argument(parser, "--threads", help="Number of Threads to Use", default=num_cores)
parser <- add_argument(parser, "--output", help="Output file", default="filtered.dataframe.tsv")
parser <- add_argument(parser, "--type", help="Whether Data is Variant/Haplotype Level", default="variant")

# Parse Arguments
args <- parse_args(parser)

# Read Input Data
ref.data <- read.table(args$refcounts, sep='\t', header=TRUE)
alt.data <- read.table(args$altcounts, sep='\t', header=TRUE)

# Convert to Numeric (To be Safe)
ref.data[, c(-2, -1)] <- sapply(ref.data[, c(-2, -1)], as.numeric)
alt.data[, c(-2, -1)] <- sapply(alt.data[, c(-2, -1)], as.numeric)


# Run QC Function
dataframe_results <- QC(ref.data, alt.data, numCores=args$threads)

# Filter good quality ASEQC data and write to File

write.table(dataframe_results, file=args$output, sep="\t")
