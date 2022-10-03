## SCRIPT NAME: 2_alevin_QC.R
## This R script generates a QC report summarizing the output of alevin (Srivastava et al., Genome Biology 20:65 (2019)).
## Package vignettte: https://bioconductor.org/packages/release/bioc/html/alevinQC.html
## NOTE: This R script needs to be sourced
## NOTE: This R script needs to be located is a sub-directory of the main directory.
## NOTE: The alevin output needs to be located in a sub-directory of the main directory called "2_salmon_alevin_output".

if(!require(rstudioapi)){
  install.packages("rstudioapi")
  library(rstudioapi)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("alevinQC")

library(alevinQC)

## Specify data and results directories
script_dir <- dirname(parent.frame(2)$ofile)  
setwd(script_dir)
maindir <- dirname(script_dir)
datadir <- paste0(maindir, "/2_salmon_alevin_output/")
resultsdir <- paste0(maindir, "/3_salmon_alevin_QC/")
dir.create(resultsdir)

alevin <- list()
for (dir in list.dirs(path = datadir, recursive = FALSE)) {
  ## Check input files
  print(checkAlevinInputFiles(dir))
  ## Run alevinQC
  alevinQCReport(baseDir = dir,
                 sampleId = basename(dir),
                 outputFile = paste0(basename(dir), "_alevinReport.html"),
                 outputDir = resultsdir, forceOverwrite = TRUE)
  ## Retrieve raw QC data
  alevin[[dir]] <- readAlevinQC(baseDir = dir)
  ## write.csv(alevin[[dir]]$cbTable, file = paste0(resultsdir, paste('cbTable', basename(dir), 'csv', sep = '.')), quote = FALSE)
  write.csv(data.frame(summaryTable = alevin[[dir]]$summaryTables$fullDataset), file = paste0(resultsdir, paste('summaryTable', basename(dir), 'csv', sep = '.')), quote = FALSE)
  write.csv(data.frame(versionTable = alevin[[dir]]$versionTable), file = paste0(resultsdir, paste('versionTable', basename(dir), 'csv', sep = '.')), quote = FALSE)
}


