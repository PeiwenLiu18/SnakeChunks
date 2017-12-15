## QUESTIONS for Claire / answers following
## - Il faut lire le fichier VolcanoPlots.R avec la commande R source()
## > À quel moment/quel endroit ? J'ai commenté le passage en attendant 
## - Tab-separated value extension: tsv or tab ? Now people use "tsv"
## > OK
## - a design file can contain several comparisons between two conditions -> one separate folder per comparison ?
## > J'ai repris la structure de dossiers que j'utilisais avec SARTools (et qui est similaire à celle du ChIP-seq)
## - design file: is there a reason for specifying the reference in the first column and the test in the second one ? This is somewhat confusing since in the output we have test_vs_ref
## > I can't remember, normally columns are read using the headers so it shouldn't matter
## - can we assume that the order of the columns in the design file is the same as rows in the sample description table ?
## > You mean the columns in the count file? I think yes but I'm not 100% sure
## - Is there a way to pass an optional list of parameters from snakemake to R in the same way as the "..." specification for function headers ?
## > I don't think I understand what you mean

## Load required librarires
library("DESeq2")


## Get current script directory
## in order to load other scripts in the same dir as this one. 
## Solution proposed Bernhard Kausler on
## https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) message("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) message("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}
script.dir <- getScriptPath()
message("Script directory: ", script.dir)
source(file.path(script.dir, "VolcanoPlot.R"))
# source("SnakeChunks/scripts/R-scripts/R/VolcanoPlot.R")


## Report working directory
message("Working directory: ", getwd())



################ Parameter values ################
## (THESE SHOULD BE PASSED VIA SNAKEMAKE)
parameters <- list()

## Get and check DESeq2 parameters

## Directory containing the differential expression
parameters$diffexpr.dir <- snakemake@params[["diffexpr_dir"]]  ## Claide, do we use this parameter ?

## Adjustment method for multiple testing
parameters$pAdjustMethod <- snakemake@params[["pAdjustMethod"]] # "BH" ## Correction for multiple testing

## Significance threshold (applied on adjusted p-value)
parameters$alpha <- as.numeric(snakemake@params[["alpha"]]) # 0.05 ## Threshold on adjusted p-value
if ((!is.numeric(parameters$alpha)) 
    || (parameters$alpha < 0) 
    || (parameters$alpha > 1)) {
  stop(parameters$alpha, " (", class(parameters$alpha), ")  is not a valid value for alpha. Must be number between 0 and 1. ")
}


## Threshold on row sums, to filter out undetected genes
parameters$rowsum_filter <- snakemake@params[["rowsum_filter"]] # 10 
if ((!is.numeric(parameters$rowsum_filter)) 
    || (parameters$rowsum_filter < 0)) {
  stop(parameters$rowsum_filter, " (", class(parameters$rowsum_filter), ")  is not a valid value for alpha. Must be a positive number. ")
}

## Data and design
parameters$sample.ids <- snakemake@params[["sample_ids"]] ## Claire: are these ever used ? 

## Tab-delimited file with sample descriptions
parameters$sample_table <- snakemake@params[["sample_tab"]] 

## Tab-delimited file with one or several pairs of conditions to be compared 
## (one comparison per row)
parameters$design_file <- snakemake@params[["design_tab"]] 

## File containing the count table (one row per feature, one column per sample)
parameters$count_file <- snakemake@input[["count_file"]] 

## Output directory
parameters$output.dir<- snakemake@params[["outdir"]]
## Create output directory if required
dir.create(parameters$output.dir, showWarnings = FALSE, recursive = TRUE)

## Read the count table
counts <- read.table(parameters$count_file, header=TRUE, row.names=1, comment.char = "#")
# dim(counts)
# View(counts)
# head(counts)

## Read sample descriptions
coldata <- read.table(parameters$sample_table, header=TRUE, row.names = 1)

## Read the design
design <- read.table(parameters$design_file, header=TRUE, row.names = NULL)


## TO CHECK: is the order of the samples the same as column names in the coutns: check with feature count rule
# colnames(counts) == rownames(coldata)
## Claire, I added a test here, please check if OK 
if (length(parameters$sample.ids) == ncol(counts)) {
  # If sample IDs are provided in user parameters, we use them
  colnames(counts) <-  parameters$sample.ids
} else if (ncol(counts) == nrow(coldata)) { 
  colnames(counts) <- rownames(coldata)
} else {
  ## Extract sample IDs from the column headers (a bit tricky)
  colnames(counts) <- sub(pattern = "RNA.seq.results.samples.", replacement = "", names(counts))
  colnames(counts) <- sub(pattern = "\\..*", replacement  = "", x = colnames(counts), perl=TRUE)
  
  ## Claire: it may be better to enforce the correct specification of sample IDs, either via the sample table or via user-specified parameters
  ## -> rather than parsing the column headers, die with an error message
  ## We should discuss this. In the meantime I put the stop. 
  stop("Invalid specification of sample IDs. ")
}

## Iterate over each line of the design file.
## Note: a design fle can contain several differential expression analyses. 
## The two first columns of each row specify the two conditions to be compared. 
i <- 1
for (i in 1:nrow(design)) {
  ## Define reference and test conditions from the design file
  ref.condition <- as.vector(unlist(design[i, 1]))
  test.condition <- as.vector(unlist(design[i, 2]))
  
  ## Select reference and test samples
  ref.samples <- row.names(coldata)[as.vector(coldata$Condition) == ref.condition]
  test.samples <- row.names(coldata)[as.vector(coldata$Condition) == test.condition]
  selected.samples <- c(ref.samples, test.samples)
  
  message ("DESEq2 analysis ", i, "/", nrow(design), 
           "; ref condition: ", ref.condition, " (", length(ref.samples)," samples)",
           "; test condition: ", test.condition, " (", length(test.samples)," samples)")  
  
  ## Build a DESeq dataset
  dds <- DESeqDataSetFromMatrix(countData=counts[, selected.samples], colData=coldata[selected.samples,], design = ~Condition)
  
  ## Define the conditions to be compared
  dds$Condition <- relevel(dds$Condition, ref=ref.condition) 
  
  # remove uninformative columns ## ???? rows, genes
  ## Filter out genes with zero counts in all samples
  dds <- dds[ rowSums(counts(dds)) > parameters$rowsum_filter, ]
  
  # Normalization, preprocessing and differential analysis
  dds <- DESeq(dds)
  
  ## Extract the result
  contrast <- c("Condition", c("FNR", "WT"))
  res <- results(dds, contrast=contrast, 
                 alpha = parameters$alpha,
                 independentFiltering=FALSE, pAdjustMethod = parameters$pAdjustMethod)  ## Collect the result table
  # dim(res)
  
  ## CLAIRE: la fonction suivante n'existe apparemment pas dans DESeq2
  # shrink fold changes for lowly expressed genes
  # res <- lfcShrink(dds, contrast=contrast, res=res) 
  # > Oui je crois qu'elle a disparu entre deux versions de DESeq2, mais ce n'est pas indispensable
  
  
  # sort by p-value
  res.sorted <- res[order(res$padj),]
  # View(data.frame(res.sorted))

  ## Build a data.frame for export
  res.frame <- cbind("gene" = row.names(res.sorted), data.frame(res.sorted))
  # names(res.frame)
  # head(res.frame)
  
  ## Select differentially expressed genes
  DEG.genes <- res.frame$padj < parameters$alpha
  
  ################ Export result files ################
  message("Exporting results to directory ", parameters$output.dir)
  ## Build prefix from conditions
  file.prefix <- file.path(parameters$output.dir, paste(sep="_", test.condition, "vs", ref.condition))
  
  ## Draw an MA plot
  pdf(snakemake@output[["ma_plot"]], width = 7, height = 7)
  DESeq2::plotMA(res.sorted, alpha=parameters$alpha, las=1)
  grid()
  silence <- dev.off()
  
  ## Volcano plot
  pdf(snakemake@output[["volcano_plot"]], width = 7, height = 7)
  VolcanoPlot(multitest.table = res.frame, 
              main=paste(sep="", test.condition, " vs ", ref.condition),
              alpha = parameters$alpha,
              effect.size.col = "log2FoldChange",
              control.type = "padj", legend.corner = "top", legend.cex = 0.8, las=1, col.positive = "#BB0000")
  silence <- dev.off()

## P-value histogram (unadjusted p-values)
pdf(snakemake@output[["pval_histo"]], width = 8, height = 6)

hist(res.frame$pvalue, main="P-value histogram", 
     breaks=seq(from=0, to=1, by=0.05),
     col="#CC6666", xlab="Nominal p-value", ylab="Number of genes", las=1)

## Highlight number of genes declared positive
hist(unlist(res.frame[!DEG.genes, "pvalue"]), 
     breaks=seq(from=0, to=1, by=0.05), add = TRUE, col="#BBBBBB")
## Estimate the number of genes under null hypothesis
n.genes <- nrow(res.frame) ## Number of genes

## Estimation of the number of genes under null (n0) or alternative (n1) hypothesis, 
## based on the method defined by Storey and Tibshirani (2003).
m <- nrow(res.frame)
n0 <- min(2*sum(res.frame$pvalue >= 0.5), m) 
n1 <- m - n0
abline(h=n0/20, lty="dashed", col="#0000BB", lwd=2)
legend("topright", 
       legend = c(
         paste("N = ", n.genes),
         paste("n0 = ", n0),
         paste("n1 = ", n1),
         paste("DEG = ", sum(DEG.genes)),
         paste("Sn = ", signif(digits=2, sum(DEG.genes) / n1))
       ), 
       lty=c("solid", "dashed", "solid", "solid", "solid"), 
       lwd=c(0,2,0,7,0), col=c(NA,"#0000BB",NA, "#CC6666", NA))
silence <- dev.off()


  ## Print a result table with all genes
  write.table(res.frame, row.names = FALSE, col.names=TRUE,
              sep="\t", quote=FALSE, 
              file=snakemake@output[["gene_table"]])
  
  ## Print a result table with genes passing the threshold
  write.table(res.frame[DEG.genes, ], row.names = FALSE, col.names=TRUE,
              sep="\t", quote=FALSE, 
#              file=paste(sep="", file.prefix, "_deseq2_DEG_", parameters$pAdjustMethod, "_alpha", parameters$alpha, ".tsv"))
              file=snakemake@output[["gene_pval"]])

  ## Export the list of differentially expressed gene names
  write.table(res.frame[DEG.genes, "gene"], row.names = FALSE, col.names=FALSE,
              sep="\t", quote=FALSE, 
#              file=paste(sep="", file.prefix, "_deseq2_DEG_", parameters$pAdjustMethod, "_alpha", parameters$alpha, "_genes.txt")) ## snakemake@output[["gene_list"]]
              file=snakemake@output[["gene_list"]])

  list.files(parameters$output.dir)
  # system(paste("open", parameters$output.dir)) ## to check the results; only works for Mac
}

