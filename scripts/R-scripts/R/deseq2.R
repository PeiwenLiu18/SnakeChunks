
#counts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
#coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample")

dir.main <- '~/FNR_analysis/'
setwd(dir.main)

library("DESeq2")

dir.counts <- 'RNA-seq/results/diffexpr'

count.file <- file.path(dir.counts, "cutadapt_bowtie2_featureCounts_all.txt")

counts <- read.table(count.file, header=TRUE, row.names="gene")

coldata <- read.table("samples_RNA-seq.tab", header=TRUE, row.names="sample")

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
