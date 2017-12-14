
#counts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
#coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample")

## QUESTIONS for Claire
## - Tab-separated value extension: tsv or tab ? Now people use "tsv"
## - a design file can contain several comparisons between two conditions -> one separate folder per comparison ?
## - design file: is there a reason for specifying the reference in the first column and the test in the second one ? This is somewhat confusing since in the output we have test_vs_ref
## - can we assume that the order of the columns in the design file is the same as rows in the sample description table ?
## - Is there a way to pass an optional list of parameters from snakemake to R in the same way as the "..." specification for function headers ?

dir.main <- '~/FNR_analysis/'
setwd(dir.main)
dir.counts <- 'RNA-seq/results/diffexpr'

################ Parameter values ################
## (THESE SHOULD BE PASSED VIA SNAKEMAKE)
parameters <- list()

## DESeq2 parameters
parameters[["pAdjustMethod"]] <- "BH" ## Correction for multiple testing
parameters[["alpha"]] <- 0.05 ## Threshold on adjusted p-value
parameters[["rowsum_filter"]] <- 10 ## Threshold on row sums, to filter out undetected genes

## Data and design
parameters[["sample.ids"]] <- "" ##
parameters[["sample_table"]] <- file.path(dir.main, "metadata", "samples_RNA-seq.tab")
parameters[["design_file"]] <- file.path(dir.main, "metadata", "design_RNA-seq.tab")
parameters[["count_file"]] <- file.path(dir.counts, "cutadapt_bowtie2_featureCounts_all.txt")

## Output directory
parameters[["output_dir"]] <- file.path(dir.counts, "DEG", "DESeq2")

## Create output directory if required
dir.output <- parameters[["output_dir"]]
dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

## The real script starts here

library("DESeq2")


counts <- read.table(parameters[["count_file"]], header=TRUE, row.names=1, comment.char = "#")
colnames(counts) <- sub(pattern = "RNA.seq.results.samples.", replacement = "", names(counts))
colnames(counts) <- sub(pattern = "\\..*", replacement  = "", x = colnames(counts), perl=TRUE)
# dim(counts)
# View(counts)
# head(counts)

## Read sample descriptions
coldata <- read.table(parameters[["sample_table"]], header=TRUE, row.names = 1)

## Read the design
design <- read.table(parameters[["design_file"]], header=TRUE, row.names = NULL)


## TO CHECK: is the order of the samples the same as column names in the coutns: check with feature count rule
# colnames(counts) == rownames(coldata)
colnames(counts) <- rownames(coldata)

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
  dds <- dds[ rowSums(counts(dds)) > parameters[["rowsum_filter"]], ]
  
  # Normalization, preprocessing and differential analysis
  dds <- DESeq(dds)
  
  ## Extract the result
  contrast <- c("Condition", c("FNR", "WT"))
  res <- results(dds, contrast=contrast, 
                 alpha = parameters[["alpha"]],
                 independentFiltering=FALSE, pAdjustMethod = "BH")  ## Collect the result table
  dim(res)
  
  ## CLAIRE: la fonction suivante n'existe apparemment pas dans DESeq2
  # shrink fold changes for lowly expressed genes
  # res <- lfcShrink(dds, contrast=contrast, res=res) 
  
  
  # sort by p-value
  res.sorted <- res[order(res$padj),]
  # View(data.frame(res.sorted))

  ## Build a data.frame for export
  res.frame <- cbind("gene" = row.names(res.sorted), data.frame(res.sorted))
  # names(res.frame)
  # head(res.frame)
  
  ## Select differentially expressed genes
  DEG.genes <- res.frame$padj < parameters[["alpha"]]
  
  ################ Export result files ################
  message("Exporting results to directory ", dir.output)
  ## Build prefix from conditions
  file.prefix <- file.path(dir.output, paste(sep="_", test.condition, "vs", ref.condition))
  
  ## Draw an MA plot
  pdf(paste(sep="_", file.prefix, "ma_plot.pdf"))
  plotMA(res.sorted)
  silence <- dev.off()
  
  
  ## Print a result table with all genes
  write.table(res.frame, row.names = FALSE, col.names=TRUE,
              sep="\t", quote=FALSE, 
              file=paste(sep="_", file.prefix, "deseq2_all_genes.tsv"))
  
  ## Print a result table with genes passing the threshold
  write.table(res.frame[DEG.genes, ], row.names = FALSE, col.names=TRUE,
              sep="\t", quote=FALSE, 
              file=paste(sep="", file.prefix, "deseq2_DEG_", parameters[["pAdjustMethod"]], "_alpha", parameters[["alpha"]], ".tsv"))
  
  ## Export the list of differentially expressed gene names
  write.table(res.frame[DEG.genes, "gene"], row.names = FALSE, col.names=FALSE,
              sep="\t", quote=FALSE, 
              file=paste(sep="", file.prefix, "deseq2_DEG_", parameters[["pAdjustMethod"]], "_alpha", parameters[["alpha"]], "_genes.txt"))
  
  list.files(dir.output)
  # system(paste("open", dir.output)) ## to check the results; only works for Mac
}

