# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
#cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
#coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample")
#print(cts)

library("DESeq2")

cts <- read.table("cutadapt_bowtie2_featureCounts_all.txt", header=TRUE, row.names="gene")
coldata <- read.table("samples_RNA-seq.tab", header=TRUE, row.names="sample")

dds <- DESeqDataSetFromMatrix(countData=cts, colData=coldata, design=~condition)

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
# TODO explore IHW usage


# store results
pdf("ma_plot.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file="deseq2_res.tab")
