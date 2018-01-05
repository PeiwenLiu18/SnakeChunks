library(VennDiagram)

chip <- snakemake@input[["chip_genes"]] #as.vector(read.table("{input.chip_genes}")[,1])
rna <- snakemake@input[["rna_genes"]] # as.vector(read.table("{input.rna_genes}")[,1])
#regulon <- as.vector(read.table("{input.regulon_genes}")[,1])

#venn.plot <- venn.diagram(list(ChIP=chip, RNA=rna, Regulon=regulon), filename="{output}", imagetype="png", fill=rainbow(3))
venn.plot <- venn.diagram(list(ChIP=chip, RNA=rna), filename="{output}", imagetype="png", fill=rainbow(3))
