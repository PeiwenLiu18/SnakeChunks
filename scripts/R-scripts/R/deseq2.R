
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
parameters[["count_file"]] <- file.path(dir.counts, "cutadapt_bowtie2_featureCounts_all.txt")


## The real script starts here

library("DESeq2")


counts <- read.table(parameters[["count_file"]], header=TRUE, row.names=1, comment.char = "#")
dim(counts)
View(counts)
coldata <- read.table(parameters[["sample_table"]], header=TRUE, row.names = 1)

# parameters[["sample.ids"]] <- sub(pattern = "RNA.seq.results.samples.", replacement = "", names(counts))
# parameters[["sample.ids"]] <- sub(pattern = "\\..*", replacement  = "", x = parameters[["sample.ids"]], perl=TRUE)

## TO CHECK: is the order of the samples the same as column names in the coutns: check with feature count rule
colnames(counts) <- rownames(coldata)


dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)

# remove uninformative columns
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
