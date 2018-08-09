
## ---- Command-line options ----

## By  default this script is invoked on the command line with options parsed from the arguments.
##
## Alternatively, the options can be defined separately, so that the script can be called from another R script.
##
## Integration in a snakemake rule will be evaluated soon.
#

if (!exists("opt")) {
  message("Reading parameters from the command line")
  ## ---- Parse command-line arguments -------
  #!/path/to/Rscript
  library('getopt')
  #get options, using the spec as defined by the enclosed list.
  #we read the options from the default: commandArgs(TRUE).
  spec = matrix(c(
    'verbose', 'v', 2, "integer",
    'help', 'h', 0, "logical",
    'main_dir', 'm', 1, "character",
    'config_file', 'c', 1, "character",
    'count_table', 't', 1, "character",
    'output_dir', 'o', 1, "character",
    'snakechunks_dir', 's', 1, 'character',
    'rsnakechunks_dir', 'r', 1, 'character'
  ), byrow = TRUE, ncol = 4)
  opt = getopt(spec)
}



## Help message
if (!is.null(opt$help) ) {
  message("\n", getopt(spec, usage = TRUE))
  q(status = 0)
}



## ---- Mandatory arguments ----
if (is.null(opt$config_file) ) {
  stop("Configuration file is a mandatory argument (-c, --config_file)")
}

if (is.null(opt$count_table) ) {
  stop("Count table is a mandatory argument (-t, --count_table)")
}

## ---- Optional arguments ----

## Verbosity
if (is.null(opt$verbose) ) {
  opt$verbose    = 1
}

## Main directory
if (is.null(opt$main_dir) ) {
  opt$main_dir = getwd()
}

## Result directory
if (is.null(opt$output_dir)) {
  opt$output_dir <- "results"
  message("Result dir defined from main directory\t", opt$output_dir)
}

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
  count.table = opt$count_table,
  configFile = opt$config_file,
  main.dir = opt$main_dir,
  result.dir = opt$output_dir,
  verbose = opt$verbose)

# script.name <- get_Rscript_filename()

# q(status = 0)
