#' @title compare gene lists returned by the ChIP-seq and RNA-seq workflow resuls
#' @author Claire Rioualen and Jacques van Helden
#' @description compare gene list returned by RNA-seq and ChIP-seq analyses:
#' \itemize{
#' \item DEG genes detected by RNA-seq
#' \item TF target genes detected by ChIp-seq
#' \item reference target genes annotated in RegulonDB
#' }

#' @examples
#' parameters <- c(
#'   chip_genes = "ChIP-seq/results/peaks/FNR1_vs_input1/macs2/FNR1_vs_input1_cutadapt_bowtie2_macs2_gene_list.txt",
#'   rna_genes = "RNA-seq/results/diffexpr/FNR_vs_WT/DESeq2/FNR_vs_WT_cutadapt_bwa_featureCounts_DESeq2_gene_list.txt",
#'   gene_descriptions = "data/RegulonDB/GeneProductSet.txt",
#'   TUs = "data/RegulonDB/TUSet.txt",
#'   BS = "data/RegulonDB/BindingSiteSet.txt",
#'   TFs = "data/RegulonDB/GeneProductSet.txt",
#'   TF = "FNR",
#'   venn.format = "png",
#'   verbose = 1,
#'   outdir = "integration"
#' )
#'
#' output <- list(
#'   "dir" = parameters[["outdir"]],
#'   "venn" = "integration/ChIP-RNA-regulons_venn.png",
#'   "gene_table" = "integration/ChIP-RNA-regulons_gene_table.tsv",
#'   "chip_genes_gff" = "integration/ChIP-RNA-regulons_ChIP-seq.gff",
#'   "rna_genes_gff" = "integration/ChIP-RNA-regulons_RNA-seq.gff",
#'   "regulon_genes_gff" = "integration/ChIP-RNA-regulons_RegulonDB.gff"
#' )
#' source("SnakeChunks/scripts/RSnakeChunks/misc/gene_lists_comparison.R")
#'

## Load libraries
library(VennDiagram)

message("Gene list comparaison: ChIP-seq vs RNA-seq vs ReegulonDB")


test.cmd <- "
Rscript SnakeChunks//scripts/RSnakeChunks/misc/gene_lists_comparison.R \
  --chip_genes=ChIP-seq/results/peaks/FNR1_vs_input1/macs2/FNR1_vs_input1_cutadapt_bowtie2_macs2_gene_list.txt \
  --rna_genes=RNA-seq/results/diffexpr/FNR_vs_WT/DESeq2/FNR_vs_WT_cutadapt_bwa_featureCounts_DESeq2_gene_list.txt \
  --gene_descriptions=data/RegulonDB/GeneProductSet.txt \
  --TUs=data/RegulonDB/TUSet.txt \
  --BS=data/RegulonDB/BindingSiteSet.txt \
  --TFs=data/RegulonDB/GeneProductSet.txt \
  --myTF=FNR \
  --venn.format=png \
  --verbose=1 \
  --outdir=integration \
  --venn=integration/ChIP-RNA-regulons_venn.png \
  --gene_table=integration/ChIP-RNA-regulons_gene_table.tsv \
  --chip_genes_gff=integration/ChIP-RNA-regulons_ChIP-seq.gff \
  --rna_genes_gff=integration/ChIP-RNA-regulons_RNA-seq.gff \
  --regulon_genes_gff=integration/ChIP-RNA-regulons_RegulonDB.gff
"



## Parameters from Snakemake rule
if (exists("parameters")) {
  message("Parameters provided within R")
} else if (exists("snakemake")) {
  message("Getting parameters from snakemake")
  parameters <- list(
    "chip_genes" = snakemake@input[["chip_genes"]],
    "rna_genes" = snakemake@input[["rna_genes"]],
    "gene_descriptions" = snakemake@input[["regulondb_gene_product"]],
    "TFBS" = snakemake@input[["regulondb_sites"]],
    "TUs" = snakemake@input[["regulondb_tus"]],
    "TF" = snakemake@params[["TF"]],
    "venn.format" = snakemake@params[["image_format"]],
    verbose = 0
  )

  output <- list(
    "dir" = snakemake@params[["outdir"]],
    "venn" = snakemake@output[["venn"]],
    "gene_table" = snakemake@output[["gene_table"]],
    "chip_genes_gff" = snakemake@output[["chip_genes_gff"]],
    "rna_genes_gff" = snakemake@output[["rna_genes_gff"]],
    "regulon_genes_gff" = snakemake@output[["regulon_genes_gff"]]
  )
} else {
  message("Reading parameters from the command line")
  ## ---- Parse command-line arguments -------
  #!/path/to/Rscript
  library('getopt')
  #get options, using the spec as defined by the enclosed list.
  #we read the options from the default: commandArgs(TRUE).
  # spec <- matrix(c(
  #   'verbose', 'v', 2, "integer",
  #   'help', 'h', 0, "logical",
  #   'chip_genes', 'c', 1, 'character',
  #   'rna_genes', 'r', 1, 'character',
  #   'gene_description', 'd', 1, 'character',
  #   'BS', 'b', 1, 'character',
  #   'TUs', 'u', 1, 'character',
  #   'TFs', 't', 1, 'character',
  #   'Venn.format', 'f', 1, 'character',
  #   "outdir", "o", 1, "character",
  #   "venn", "e", 1, "character",
  #   "gene_table", "a", 1, "character",
  #   "chip_genes_gff", "d", 1, "character",
  #   "rna_genes_gff", "e",1, "character",
  #   "regulon_genes_gff", "i",1, "character"
  #   # 'snakechunks_dir', 's', 1, 'character',
  #   # 'rsnakechunks_dir', 'r', 1, 'character',
  # ), byrow = TRUE, ncol = 4)
  # opt = getopt(spec)


  library("argparse")
  parser <- ArgumentParser()
#  parser$add_argument("-h", "--help", action=store_false,  degault = TRUE, help = "Get help message. ")
  parser$add_argument("-v", "--verbose", default = 0, type = "integer", help = "verbosity level")
  parser$add_argument("-c", "--chip_genes", default = "", type = "character", help = "File containing ChIP-seq peak-associated genes.")
  parser$add_argument("-r", "--rna_genes", default = "", type = "character", help = "File containing differentially expressed genes from RNA-seq.")
  parser$add_argument("-d", "--gene_descriptions", default = "", type = "character", help = "File containing a table with gene descriptions.")
  parser$add_argument("-b", "--BS", default = "", type = "character", help = "Table of reference transcriptiojn factor binding sites (from RegulonDB).")
  parser$add_argument("-u", "--TUs", default = "", type = "character", help = "Table of transcription units (from RegulonDB).")
  parser$add_argument("-t", "--TFs", default = "", type = "character", help = "Table of transcription factors (from RegulonDB).")
  parser$add_argument("-m", "--myTF", default = "", type = "character", help = "Transcription factor of interest.")
  parser$add_argument("-o", "--outdir", default = "", type = "character", help = "output directory.")
  parser$add_argument("-n", "--venn", default = "", type = "character", help = "output file with the Venn diagram")
  parser$add_argument("-e", "--venn.format", default = "", type = "character", help = "output file format for the Venn diagram (supported: png, pdf).")
  parser$add_argument("-a", "--gene_table", default = "", type = "character", help = "output file with the annotated genes.")
  parser$add_argument("-g", "--chip_genes_gff", default = "", type = "character", help = "output GFF-formatted file with the ChIP-seq peak associated genes (for IGV viewer).")
  parser$add_argument("-x", "--rna_genes_gff", default = "", type = "character", help = "output GFF-formatted file with the RNA-seq DEG (for IGV viewer).")
  parser$add_argument("-l", "--regulon_genes_gff", default = "", type = "character", help = "output GFF-formatted file with the know regulon (for IGV viewer).")
  opt <- parser$parse_args()

  # ## Help message
  # if (!is.null(opt$help) ) {
  #   message("\n", getopt(spec, usage = TRUE))
  #   q(status = 0)
  # }

  parameters <- list(
    "chip_genes" = opt[["chip_genes"]],
    "rna_genes" = opt[["rna_genes"]],
    "gene_descriptions" = opt[["gene_descriptions"]],
    "BS" = opt[["BS"]],
    "TFs" = opt[["TFs"]],
    "TUs" = opt[["TUs"]],
    "TF" = opt[["myTF"]],
    "venn.format" = opt[["venn.format"]],
    verbose = 0
  )
  #message("Parameters\t", (parameters))

  output <- list(
    "dir" = opt[["outdir"]],
    "venn" = opt[["venn"]],
    "gene_table" = opt[["gene_table"]],
    "chip_genes_gff" = opt[["chip_genes_gff"]],
    "rna_genes_gff" = opt[["rna_genes_gff"]],
    "regulon_genes_gff" = opt[["regulon_genes_gff"]]
  )

  #message("Output optiont\t", t(as.data.frame(output)))

  ## ---- Mandatory arguments ----
  if (is.null(opt$chip_genes) ) {
    stop("ChIP-seq gene list is a mandatory argument (-c, --chip_genes)")
  }
  parameters[["chip_genes"]] <- opt$chip_genes

  if (is.null(opt$rna_genes) ) {
    stop("RNA-seq gene list is a mandatory argument (-r, --rna_genes)")
  }
  parameters[["rna_genes"]] <- opt$rna_genes

  if (is.null(opt$gene_description) ) {
    stop("Gene description table is a mandatory argument (-g, --gene_description)")
  }
  parameters[["gene_descriptions"]] <- opt$gene_descriptions


}


## ---- Optional arguments ----
if (is.null(parameters[["verbose"]])) {
  parameters[["verbose"]] = 1
}

## Report parameters
if (parameters[["verbose"]] >= 0) {
  message("\tParameters")
  for (p in names(parameters)) {
    message("\t\t", p, "\t", parameters[p])
  }

  message("\tOutput options")
  for (p in names(output)) {
    message("\t\t", p, "\t", output[p])
  }
}



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
message("\tReading gene description table\t", parameters[["gene_descriptions"]])
gene.table <- read.delim(file = parameters[["gene_descriptions"]],
                         comment.char = "#", as.is = TRUE,
                         quote = NULL)
names(gene.table) <- c("gene_id", "gene_name", "bnumber", "gene_left", "gene_right", "strand",
                       "product", "evidence", "PIMDs", "evidence_level")
message("\t\t", nrow(gene.table), " genes")
# View(gene.table)

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
message("\tReading transcription factor binding sites (TFBS) from file\t", parameters[["BS"]])
TFBS <- read.delim(file = parameters[["BS"]],
                   comment.char = "#",
                   quote = NULL)
names(TFBS) <- c("TF_id", "TF_name", "TFBS_id", "TFBS_left", "TFBS_right", "strand",
                 "interaction_id", "TU_name", "effect", "promoter_name", "TFBS_center", "TFBS_sequence",
                 "evidence", "conf_level")
message("\t\t", nrow(TFBS), " TFBS")
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
message("\tReading transcription units (TUs) from file\t", parameters[["TUs"]])
TUs <- read.delim(file = parameters[["TUs"]],
                  comment.char = "#",
                  quote = NULL)
message("\t\t", nrow(TUs), " TUs")
names(TUs) <- c("TU_id", "TU_name", "operon_name", "gene_names", "promoter_name", "evidence", "conf_level")
# View(TUs)



## ---- Get info for TF of interest ----
message("\tGetting RegulonDB data for factor ", parameters[["TF"]])

#### Select reference sites ####
ref.sites <- subset(TFBS, TF_name == parameters[["TF"]])
if (nrow(ref.sites) == 0) {
  stop("RegulonDB does not contain any binding site for transcription factor ", parameters[["TF"]])
}
message("\t\t", nrow(ref.sites), " TFBS")

#### Identify target transcription units ####
target.TUs <- unique(sort(as.vector(unlist(ref.sites$TU_name))))
# names(ref.sites)
message("\t\t", length(target.TUs), " TUs")
# message(paste(collapse = ", ", target.TUs))

#### Get target genes from the TFBS ####
tu.gene.names <- unique(sort(as.vector(as.matrix(subset(TUs, TU_name %in% target.TUs, select = "gene_names")))))
message("\t\t", length(tu.gene.names), " target TU gene names\t")
# message(paste(collapse = ", ", tu.gene.names))

target.genes  <-  unique(sort(unlist(strsplit(x = tu.gene.names, split = ",", fixed = TRUE))))
message("\t\t", length(target.genes), " target genes")
target.gene.ids <- unlist(subset(gene.table, gene_name %in% target.genes, select = bnumber))
message("\t\t", length(target.gene.ids), " target gene IDs")

#### Load gene lists ####
genes <- list(
  "ChIPseq" = scan(parameters[["chip_genes"]], what = "character", quiet = TRUE),
  "RNAseq" =  scan(parameters[["rna_genes"]], what = "character", quiet = TRUE),
  "RegulonDB" = target.gene.ids
)


## Create output directory
dir.create(output[["dir"]], showWarnings = FALSE, recursive = TRUE)

#venn.plot <- venn.diagram(list(ChIP=chip, RNA=rna, Regulon=regulon), filename="{output}", imagetype="png", fill=rainbow(3))
#for (venn.format in unlist(strsplit(parameters[["venn.format"]], split = " "))) {
venn.file <- output[["venn"]]
message("Exporting Venn diagram ", venn.file)
venn.plot <- venn.diagram(genes,
                          filename = venn.file,
                          imagetype = parameters[["venn.format"]],
                          fill = rainbow(length(genes)))

#}

#### Export summary table with the different criteria ####
row.names(gene.table) <- gene.table$gene_id
regulon.name <- paste(parameters[["TF"]], "regulon", sep = "_")
gene.table[, "ChIPseq"] <- 0
gene.table[, "RNAseq"] <- 0
gene.table[, regulon.name] <- 0
gene.table[gene.table$bnumber %in% genes$ChIPseq, "ChIPseq"] <- 1
gene.table[gene.table$bnumber %in% genes$RNAseq, "RNAseq"] <- 1
gene.table[gene.table$bnumber %in% genes$RegulonDB, regulon.name] <- 1
out.gene.table <- output[["gene_table"]]
message("Exporting annotated gene table: ", out.gene.table)
write.table(x = gene.table, sep = "\t", quote = FALSE,
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
  "seqname" = "Chromosome",
  "source" = "SnakeChunks",
  "feature" = "gene",
  "start" = gene.table$gene_left,
  "end" = gene.table$gene_right,
  "score" = ".",
  "strand" = sub(pattern ="reverse", replacement = "-", sub(pattern = "forward", replacement = "+", x = gene.table$strand)),
  "frame" = ".",
  "attribute" = paste(sep = "", "gene_id: ", gene.table$bnumber)
)

chipseq.gff <- output[["chip_genes_gff"]]
message('Exporting GFF file for ChIP-seq results: ', chipseq.gff)
write.table(x = subset(gff, gene.table$ChIPseq == 1), file = chipseq.gff, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

rnaseq.gff <- output[["rna_genes_gff"]]
message('Exporting GFF file for RNA-seq results: ', rnaseq.gff)
write.table(x = subset(gff, gene.table$RNAseq == 1), file = rnaseq.gff, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

regulon.gff <- output[["regulon_genes_gff"]]
message('Exporting GFF file for RegulonDB results: ', regulon.gff)
write.table(x = subset(gff, gene.table[[regulon.name]] == 1), file = regulon.gff, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

