## compare gene list returned by different analyses
## - DEG genes detected by RNA-seq
## - TF target genes detected by ChIp-seq
## - reference target genes annotated in RegulonDB

library(VennDiagram)

dir.main <- "~/FNR_analysis"
setwd(dir.main)

parameters <- list(
  "chip_genes" = "ChIP-seq/results/peaks/FNR_vs_input/homer/FNR_vs_input_cutadapt_bowtie2_homer_gene_list.txt",
  "rna_genes" = "RNA-seq/results/diffexpr/FNR_vs_WT/DESeq2/FNR_vs_WT_cutadapt_bwa_featureCounts_DESeq2_gene_list.txt",
  "gene.descriptions" = "data/regulonDB/GeneProductSet.txt",
  "TFBS" = "data/regulonDB/BindingSiteSet.txt",
  "TUs" = "data/regulonDB/TUSet.txt",
  "TF" = "FNR",
  "venn.format" = "png svg"
)

output <- list(
  "dir" = "integration",
  "venn" = "ChIP-RNA-regulons_Venn",
  "annotated_genes" = "ChIP-RNA-regulons_table"
)

#### Load gene description table ####
# Columns:
# (1) Gene identifier assigned by RegulonDB
# (2) Gene name
# (3) Blattner number (bnumber) of the gene
# (4) Gene left end position in the genome
# (5) Gene right end position in the genome
# (6) DNA strand where the gene is coded
# (7) Product name of the gene
# (8) Evidence that supports the existence of the gene
# (9) PMIDs list
# (10) Evidence confidence level (Confirmed, Strong, Weak)
gene.table <- read.delim(file = parameters[["gene.descriptions"]], 
                         comment.char = "#", as.is=TRUE,
                         quote = NULL)
names(gene.table) <- c("gene_id", "gene_name", "bnumber", "gene_left", "gene_right", "strand", 
                       "product", "evidence", "PIMDs", "evidence_level")
# View(gene.table)

#### Load TFBS ####
# Columns:
# (1) Transcription Factor (TF) identifier assigned by RegulonDB
# (2) TF name
# (3) TF binding site (TF-bs) identifier assigned by RegulonDB 
# (4) TF-bs left end position in the genome 
# (5) TF-bs right end position in the genome
# (6) DNA strand where the  TF-bs is located
# (7) TF-Gene interaction identifier assigned by RegulonDB (related to the "TF gene interactions" file) 
# (8) Transcription unit regulated by the TF
# (9) Gene expression effect caused by the TF bound to the  TF-bs (+ activation, - repression, +- dual, ? unknown)
# (10) Promoter name
# (11) Center position of TF-bs, relative to Transcription Start Site
# (12) TF-bs sequence (upper case)
# (13) Evidence that supports the existence of the TF-bs
# (14) Evidence confidence level (Confirmed, Strong, Weak)
TFBS <- read.delim(file = parameters[["TFBS"]], 
                         comment.char = "#",
                         quote = NULL)
names(TFBS) <- c("TF_id", "TF_name", "TFBS_id", "TFBS_left", "TFBS_right", "strand", 
                 "interaction_id", "TU_name", "effect", "promoter_name", "TFBS_center", "TFBS_sequence", 
                 "evidence", "conf_level")
# View(TFBS)


#### Load transcription units ####
# Columns:
# (1) Transcription Unit identifier assigned by RegulonDB
# (2) Transcription unit name 
# (3) Operon name containing the transcription unit
# (4) Name of the gene(s) contained in the transcription unit
# (5) Promoter Name
# (6) Evidence that supports the existence of the transcription unit
# (7) Evidence confidence level (Confirmed, Strong, Weak)
TUs <- read.delim(file = parameters[["TUs"]], 
                  comment.char = "#",
                  quote = NULL)
names(TUs) <- c("TU_id", "TU_name", "operon_name", "gene_names", "promoter_name", "evidence", "conf_level")
# View(TUs)


message("Getting RegulonDB data for factor ", parameters[["TF"]])

#### Select reference sites ####
ref.sites <- subset(TFBS, TF_name == parameters[["TF"]])
if (nrow(ref.sites) == 0) {
  stop("RegulonDB does not contain any binding site for transcription factor ", parameters[["TF"]])
}
message("\t", nrow(ref.sites), " TFBS")

#### Identify target transcription units ####
target.TUs <- unique(sort(as.vector(ref.sites$TU_name)))
message("\t", length(target.TUs), " TUs")

#### Get target genes from the TFBS ####
target.genes  <-  unique(sort(unlist(strsplit(x = as.vector(unlist(subset(TUs, TU_name %in% target.TUs, select = "gene_names"))), split = ",", fixed=TRUE))))
message("\t", length(target.genes), " target genes")

target.gene.ids <- unlist(subset(gene.table, gene_name %in% target.genes, select=bnumber))

#### Load gene lists ####
genes <- list(
  "ChIPseq" = scan(parameters[["chip_genes"]], what = "character"),
  "RNAseq" =  scan(parameters[["rna_genes"]], what = "character"),
  "regulon" = target.gene.ids
)


## Create output directory
dir.create(output[["dir"]], showWarnings = FALSE, recursive = TRUE)

#venn.plot <- venn.diagram(list(ChIP=chip, RNA=rna, Regulon=regulon), filename="{output}", imagetype="png", fill=rainbow(3))
for (venn.format in unlist(strsplit(parameters[["venn.format"]], split = " "))) {
  venn.file <- file.path(output[["dir"]], paste(sep=".", output[["venn"]], venn.format))
  message("Exporting Venn diagram", venn.file)
  venn.plot <- venn.diagram(genes, 
                            filename = venn.file, 
                            imagetype = venn.format,
                            fill=rainbow(length(genes)))

}

#### Export summary table with the different criteria ####
row.names(gene.table) <- gene.table$gene_id
regulon.name <- paste(parameters[["TF"]], "regulon", sep="_")
gene.table[, "ChIPseq"] <- 0
gene.table[, "RNAseq"] <- 0
gene.table[, regulon.name] <- 0
gene.table[gene.table$bnumber %in% genes$ChIPseq, "ChIPseq"] <- 1
gene.table[gene.table$bnumber %in% genes$RNAseq, "RNAseq"] <- 1
gene.table[gene.table$bnumber %in% genes$regulon, regulon.name] <- 1
out.gene.table <- file.path(
  output[["dir"]], 
  paste(sep=".", output[["annotated_genes"]], "tsv"))
message("Exporting annotated gene table: ", out.gene.table)
write.table(x = gene.table, sep="\t", quote=FALSE,
            row.names = FALSE,
            file = out.gene.table)


##### Export gff #####
##
## Format specifications: https://genome.ucsc.edu/FAQ/FAQformat.html#format3
## seqname - The name of the sequence. Must be a chromosome or scaffold.
## source - The program that generated this feature.
## feature - The name of this type of feature. Some examples of standard feature types are "CDS" "start_codon" "stop_codon" and "exon"li>
##   start - The starting position of the feature in the sequence. The first base is numbered 1.
## end - The ending position of the feature (inclusive).
## score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). If there is no score value, enter ":.":.
## strand - Valid entries include "+", "-", or "." (for don't know/don't care).
## frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be ".".
## group - All lines with the same group are linked together into a single item.

gff <- data.frame(
  "seqname" = "Chromsome",
  "source" = "SnakeChunks",
  "feature" = "gene",
  "start" = gene.table$gene_left,
  "end" = gene.table$gene_right,
  "score" = ".",
  "strand" = sub(pattern="reverse", replacement = "-", sub(pattern = "forward", replacement = "+", x = gene.table$strand)),
  "frame" = ".",
  "attribute" = paste(sep="", "gene_id: ", gene.table$bnumber)
)

chipseq.gff <- file.path(output[["dir"]], paste(sep="", output[["annotated_genes"]], "_ChIP-seq.tsv"))
message('Exporting GFF file for ChIP-seq results: ', chipseq.gff)
write.table(x = subset(gff, gene.table$ChIPseq == 1), file = chipseq.gff, row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)

rnapseq.gff <- file.path(output[["dir"]], paste(sep="", output[["annotated_genes"]], "_RNA-seq.tsv"))
message('Exporting GFF file for RNA-seq results: ', rnaseq.gff)
write.table(x = subset(gff, gene.table$RNAseq == 1), file = rnaseq.gff, row.names = FALSE, col.names = FALSE, sep="\t", quote=FALSE)
