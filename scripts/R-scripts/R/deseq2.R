
#counts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
#coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample")

## QUESTIONS for Claire
## - Tab-separated value extension: tsv or tab ? Now people use "tsv"
## - a design file can contain several comparisons between two conditions -> one separate folder per comparison ?

dir.main <- '~/FNR_analysis/'
setwd(dir.main)
dir.counts <- 'RNA-seq/results/diffexpr'

## Parameter values
parameters <- list()
parameters[["pAdjustMethod"]] <- "BH" ## Correction for multiple testing
parameters[["alpha"]] <- 0.05 ## Threshold on adjusted p-value
parameters[["sample.ids"]] <- ""
parameters[["sample_table"]] <- file.path(dir.main, "metadata", "samples_RNA-seq.tab")
parameters[["design_file"]] <- file.path(dir.main, "metadata", "design_RNA-seq.tab")
parameters[["count_file"]] <- file.path(dir.counts, "cutadapt_bowtie2_featureCounts_all.txt")
parameters[["dir_output"]] <- file.path(dir.counts, "DEG", "DESeq2")

## Create output directory if required
dir.output <- parameters[["dir_output"]]
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

## Iterate over each line of the design file
i <- 1
for (i in 1:nrow(design)) {
  ## Select samples according to the design
  ref.condition <- as.vector(unlist(design[i, 1]))
  test.condition <- as.vector(unlist(design[i, 2]))
  
  ref.samples <- row.names(coldata)[as.vector(coldata$Condition) == ref.condition]
  test.samples <- row.names(coldata)[as.vector(coldata$Condition) == test.condition]
  
  message ("Analysis ", i, "/", nrow(design), 
           "; ref condition: ", ref.condition, " (", length(ref.samples)," samples)",
           "; test condition: ", test.condition, " (", length(test.samples)," samples)")  
  
  selected.samples <- c(ref.samples, test.samples)
  dds <- DESeqDataSetFromMatrix(countData=counts[, selected.samples], colData=coldata[selected.samples,], design = ~Condition)
  dds$Condition <- relevel(dds$Condition, ref=ref.condition) 
  
  # remove uninformative columns ## ???? rows, genes
  ## Filter out genes with zero counts in all samples
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  
  # normalization and preprocessing
  dds <- DESeq(dds)
  
  contrast <- c("Condition", c("FNR", "WT"))
  res <- results(dds, contrast=contrast, 
                 alpha = parameters[["alpha"]],
                 independentFiltering=FALSE, pAdjustMethod = "BH")  ## Collect the result table
  dim(res)
  
  ## CLAIRE: la fonction suivante n'existe apparemment pas dans DESeq2
  # shrink fold changes for lowly expressed genes
  # res <- lfcShrink(dds, contrast=contrast, res=res) 
  
  
  # sort by p-value
  res <- res[order(res$padj),]
  View(data.frame(res))
  
  ## Build prefix from conditions
  file.prefix <- file.path(dir.output, paste(sep="_", test.condition, "vs", ref.condition))
  
  # MA plot
  pdf(paste(sep="_", file.prefix, "ma_plot.pdf"))
  plotMA(res)
  silence <- dev.off()
  
  ## Build a data.frame for export
  res.frame <- cbind("gene" = row.names(res), data.frame(res))
  # names(res.frame)
  # head(res.frame)
  
  ## Print a result table with all genes
  write.table(res.frame, row.names = FALSE, col.names=TRUE,
              sep="\t", quote=FALSE, 
              file=paste(sep="_", file.prefix, "deseq2_all_genes.tsv"))
  
  ## Print a result table with genes passing the threshold
  write.table(res.frame, row.names = FALSE, col.names=TRUE,
              sep="\t", quote=FALSE, 
              file=paste(sep="_", file.prefix, "deseq2_all_genes.tsv"))
  
  list.files(dir.output)
  # system(paste("open", dir.output)) ## to check the results; only works for Mac
}

