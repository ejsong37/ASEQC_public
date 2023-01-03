# ASEQC

**A**llele **S**pecific **E**xpression **Q**uality **C**ontrol.

This script has been authored by Eric Song and Kaushik Ganapathy from the PejLab at Scripps Research, La Jolla, California. See our paper for methods and descriptions. 


### Installation Instructions
***

* The following are dependencies of the package. To install, run the following R code in an interactive Rstudio session or terminal
``` R
install.packages(c("pracma", "data.table", "boot", "parallel", "dplyr", 
"tictoc", "stringr", "optimParallel", "R.utils", "microbenchmark", 
"argparser", "robustbase", "foreach", "doParallel", "progress"))
```

* Clone this [repository](https://github.com/Krganapa/ASEQC/) and navigate to the ASEQC folder.  

### Running Instructions
***

Run the ASEQC Script using the command
``` R
Rscript aseqc.R --refcounts <path to tsv file of ref counts> --altcounts <path to tsv file of alt counts> --threads <insert number of threads to use> --output <insert output file>
```

### Input Files
***

Both refcounts and altcounts are tab-separated (TSV) files with an index column and a column `name` containing the names of genes in addition to the columns for each sample. Ref counts will have the reference allele counts, and Alt counts will have the alternative allele counts from ASE Data. Use the included prepare_ase_input.py file to create compatible input files if necessary 

The format is outlined below:

|   | name   | sample 1 | sample 2 | ... | sample n |
| - | ------ | -------- | -------- | --- | -------- |
| 0 | Gene 1 | Count 1  | Count 2  | ... | Count n  |
| 1 | Gene 2 | Count 3  | Count 4  | ... | Count m  |
| ... | ... | ... | ... | ... | ... |

Both the index (unnamed) column and name column are **essential** for performing ASEQC fits.



