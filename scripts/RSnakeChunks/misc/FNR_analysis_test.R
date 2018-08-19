## Test to run the RNA-seq_analysis.R script in interactive R interface, with the FNR analysis study case
##
##
## Equivalent ommand-line
## Rscript --vanilla  SnakeChunks/scripts/RSnakeChunks/misc/rna-seq_deg.R -v 1 \
##    --config_file=metadata/config_RNA-seq.yml  \
##    --count_table=RNA-seq/results/diffexpr/cutadapt_bwa_featureCounts_all.tsv  \
##    --output_dir=RNA-seq/results/diffexpr/  \
##    --rsnakechunks_dir=SnakeChunks/scripts/RSnakeChunks  \
##    --report RNA-seq/results/diffexpr/cutadapt_bwa_featureCounts_rna-seq_deg_report.Rmd  \
##    &> RNA-seq/results/diffexpr/cutadapt_bwa_featureCounts_rna-seq_deg.log

analysis.dir <- '~/FNR_analysis/'

setwd(analysis.dir)

opt <- list(
  main_dir = analysis.dir,
  config_file = "metadata/config_RNA-seq.yml",
  count_table = "RNA-seq/results/diffexpr/cutadapt_bwa_featureCounts_all.tsv",
  rsnakechunks_dir = "SnakeChunks/scripts/RSnakeChunks",
  output_dir  = "RNA-seq/results/diffexpr/",
  rmd_report  = "RNA-seq/results/diffexpr/cutadapt_bwa_featureCounts_rna-seq_deg_report.Rmd",
  prefix = "cutadapt_bwa_featureCounts_",
  verbose = 1
)



## ----  Load R functions from SnakeChunks directory if specified ----
if ((!is.null(opt$snakechunks_dir)) & (is.null(opt$rsnakechunks_dir))) {
  opt$rsnakechunks_dir <- file.path(opt$snakechunks_dir, "scripts", "RSnakeChunks")
}
if (!is.null(opt$rsnakechunks_dir)) {
  message("Reading R functions from RSnakeChunks\t", opt$SnakeChunks)

  ## NOTE: if the RSnakechunks package has been compiled on this machine
  ## the functions will be loaded with library('RSnakeChunks'). However
  ## for the ime being we do not assume that RSnakeChunks has been
  ## compiled with SnakeChunks installation, so we souce all the files
  ## containing R functions.
  R.dir <- file.path(opt$rsnakechunks, "R")
  message("R directory\t", R.dir)
  R.files <- list.files(R.dir)
  for (f in R.files) {
    message("\tLoading R file ", f)
    source(file.path(R.dir, f))
  }
}

## ---- Run the analysis ----
RNAseqAnalysis(
  countFile = opt$count_table,
  configFile = opt$config_file,
  main.dir = opt$main_dir,
  result.dir = opt$output_dir,
  rmd.report = opt$rmd_report,
  prefix = opt$prefix,
  verbose = opt$verbose)

