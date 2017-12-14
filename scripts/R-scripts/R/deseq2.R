
#counts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
#coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample")


dir.main <- '~/FNR_analysis/'
setwd(dir.main)
dir.counts <- 'RNA-seq/results/diffexpr'

## Parameter values
parameters <- list()
parameters[["alpha"]] <- 0.05
parameters[["sample.ids"]] <- ""
parameters[["sample_table"]] <- file.path(dir.main, "metadata", "samples_RNA-seq.tab")
parameters[["design_file"]] <- file.path(dir.main, "metadata", "design_RNA-seq.tab")
parameters[["count_file"]] <- file.path(dir.counts, "cutadapt_bowtie2_featureCounts_all.txt")


################################################################
## DESeq2 analysis
## 
## Detect differentially expressed genes (DEG) using the package DESeq2, 
## and add a few custom columns (e-value, ...) to the result table. 
deseq2.analysis <- function(dir.figures=NULL) {
  verbose("\t\tDESeq2 analysis", 2)
  
  ## Path prefix to save DESeq2 result files
  prefix["DESeq2_file"] <- paste(sep="", prefix["comparison_file"], "_", suffix.DESeq2)
  prefix["DESeq2_figure"] <- paste(sep="", prefix["comparison_figure"], "_", suffix.DESeq2)
  
  ## Create a DESeqDataSet object from the count table + conditions
  condition <- as.factor(as.vector(current.sample.conditions))
  deseq2.dds <- DESeqDataSetFromMatrix(
    countData = current.counts, 
    colData = DataFrame(condition),
    ~ condition)
  
  
  ## Indicate that second condition is the reference condition. 
  ## If not done, the conditions are considered by alphabetical order, 
  ## which may be misleading to interpret the log2 fold changes. 
  deseq2.dds$condition <- relevel(deseq2.dds$condition, ref=cond2) 
  
  ## Run the differential analysis
  deseq2.dds <- DESeq(deseq2.dds)      ## Differential analysis with negbin distrib
  deseq2.res <- results(deseq2.dds, independentFiltering=FALSE, pAdjustMethod = "BH")  ## Collect the result table
  
  deseq2.result.table <- data.frame(
    "gene.id" = row.names(deseq2.res),
    "mean" = deseq2.res$baseMean,
    "log2FC" = deseq2.res$log2FoldChange,
    "pvalue" = deseq2.res$pvalue,
    "padj" = deseq2.res$padj)
  deseq2.result.table <- complete.deg.table(
    deseq2.result.table, 
    paste(sep="_", "DESeq2", prefix["comparison"]),
    sort.column = "padj",
    thresholds=thresholds,
    round.digits = 3,
    dir.figures=dir.figures)
  return(deseq2.result.table)
}  


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


dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design = ~Condition)

# remove uninformative columns ## ???? rows, genes
## Filter out genes with zero counts in all samples
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds)

contrast <- c("condition", c("FNR", "WT"))
res <- results(dds, contrast=contrast)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res) 
# sort by p-value
res <- res[order(res$padj),]



# store results
pdf("ma_plot.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file="deseq2_res.tab")
