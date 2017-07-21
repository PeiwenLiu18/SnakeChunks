SnakeChunks library
================================================================

The library contains a variety of files, including scripts and 
configuration files. 

You can find a description of these hereafter.

Snakemake files (snakefiles)
----------------------------------------------------------------

Snakefiles are based on the scripting language Python 3, and use a specific syntax.

For organization purpose, we have distinguished two types 
of Snakefiles: 

* Rules are typically "bricks" to build workflows with. Each rule corresponds to a specific operation.

* Workflows are combinations of rules that serve a specific purpose: quality check of sequencing data, ChIP-seq peaks analysis...


Workflows (.wf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

File extension: *.wf

import_from_sra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

quality_control
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ChIP-seq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

RNA-seq_DEG
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

integration_ChIP_RNA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Rules (.rules)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section describes briefly each rule and its parameters. 
For details on the tools and how to install them, please check `this section <http://snakechunks.readthedocs.io/en/latest/dependencies.html#>`__.

annotate_peaks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This rule runs a program that is part of the HOMER tools suite. 
It outputs a list of gene identifiers using a bed file, a fasta file and a gtf file. 
The bed file can be a peak file produces by any peak-calling rule. 

More: http://homer.salk.edu/homer/ngs/annotation.html

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["fasta_file"]
- config["genome"]["gtf_file"]

bam_by_name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sort aligned reads (in bam format) by name using 'samtools sort'. 

Required parameters:

- config["qsub"]

bam_by_pos
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sort aligned reads (in bam format) by positions, using 'samtools sort'. 

Required parameters:

- config["qsub"]

bam_stats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Computes mapping statistics using the 'samtools flagstat' tool.

Requires samtools 1.3+ version (not in apt-get repository as of 2016-03).

Required parameters:

- config["qsub"]


bam_to_bed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Converts bam files into bed files using 'bedtools bamtobed'.

Required parameters:

- config["qsub"]


bbduk
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Performs trimming of raw reads using bbduk of the bbmap suite. 
Currently only handling single-end data. 

Required parameters:

- config["qsub"]
- config["metadata"]["seq_type"]

Optional parameters:

- config["bbduk"]["length_threshold"]
- config["bbduk"]["qual_threshold"]
- config["metadata"]["strands"]


bed_to_fasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get a fasta file from a bedfile using the 'fetch-sequences' tool from the RSAT suite. 

Fetch sequences from UCSC to obtain a fasta file from a set of 
genomic coordinates described in a bed file. Requires RSAT installation. 

Alternative is to use getfasta.rules .

Example: 

::

    mkdir -p test/fetch_seq; cd test/fetch_seq; 
    wget http://pedagogix-tagc.univ-mrs.fr/rsat/demo_files/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.bed; 
    cd -
    snakemake --snakefile ${RSAT}/snakemake_files/chip-seq_motifs.py test/fetch_seq/fetch-sequences_Schmidt_2011_mm9_CEBPA_SWEMBL_R0.12_702peaks.fasta

Required parameters:

- config["qsub"]

bedgraph_to_bigwig
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Convert bedgraph to bigWig format using Deeptools. 

Required parameters:

- config["qsub"]

bedgraph_to_tdf
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Convert bedgraph to TDF format, which is recommended to load coverage data in IGV.

The conversion relies on `igvtools <https://www.broadinstitute.org/software/igv/igvtools>`__.

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["fasta_file"]


bedtools_closest
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bedtools closest searches for overlapping features in two coordinate files. 
In the event that no feature in B overlaps the current feature in A, closest will report the nearest 
(that is, least genomic distance from the start or end of A) feature in B. 

Usage: bedtools closest [OPTIONS] -a <FILE> -b <FILE1, FILE2, ..., FILEN>

More: http://bedtools.readthedocs.io/en/latest/content/tools/closest.html

Required parameters: 

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["gff3_file"]

bedtools_intersect
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bedtools intersect allows one to screen for overlaps between two sets of genomic features. 
Moreover, it allows one to have fine control as to how the intersections are reported. 
bedtools intersect works with both BED/GFF/VCF and BAM files as input.

More: http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html

Required parameters: 

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["gff3_file"]

bedtools_window
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Similar to bedtools intersect, window searches for overlapping features in A and B. 
However, window adds a specified number (1000, by default) of base pairs upstream 
and downstream of each feature in A. In effect, this allows features in B that are 
near features in A to be detected.

More: http://bedtools.readthedocs.io/en/latest/content/tools/window.html

Required parameters: 

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["gff3_file"]

Opional parameters:

- config["bedtools"]["window"]

blast_formatdb
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Run the formatdb program of the BLAST suite in order to index all 
k-mers of the reference database. This has to be done only once, then 
the DB can be used for multiple searches with blastall.

Required parameters:

- config["qsub"]

blastall
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Where {blast_program} should be replaced by one of the supported 
program options in blastall: blastp, blastn, blastx, tblastn.

Output file name: {query}_{blast_program}_hits.txt

Required parameters:

- config["qsub"]
- config["blastall"]["db"]

Optional parameters:

- config["blastall"]["matrix"]
- config["blastall"]["expect"]
- config["blastall"]["view"]


bowtie_index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rule for the creation of Bowtie 1 index. Has to be done only once.
The output file is used to test whether the index already exists when aligning.

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["version"]
- config["genome"]["fasta_file"]

bowtie
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read mapping using bowtie. 
Requires the indexing to have previously been done (using the rule bowtie_index).

Required parameters:

- config["genome"]["version"]
- config["genome"]["fasta_file"]
- config["qsub"]
- config["dir"]["fastq"]
- config["dir"]["samples"]

Optional parameters:

- config["bowtie"]["max_mismatches"]
- config["bowtie"]["threads"]


bowtie2_index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rule for the creation of Bowtie 2 index. Has to be done only once.
The output file is used to test whether the index already exists when aligning.

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["fasta_file"]

bowtie2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read mapping using Bowtie2. 
Requires the indexing to have previously been done (using the rule bowtie2_index).

Required parameters:

- config["genome"]["version"]
- config["genome"]["fasta_file"]
- config["qsub"]
- config["dir"]["fastq"]
- config["dir"]["samples"]

Optional parameters:

- config["bowtie2"]["threads"]
- config["bowtie2"]["max_mismatches"]


bPeaks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

bPeaks is a peak-calling tool running in R.
It was specifically designed for small eucaryotic organisms, such as the yeast.
It is thus not recommanded for bigger genomes, as it could be very slow. 
You should choose parameters carefully. Input in bam, output in bed.

Required parameters:

- config["qsub"]
- config["dir"]["samples"]
- config["dir"]["peaks"]

Optional parameters:

- config["bPeaks"]["IPcoeff"]
- config["bPeaks"]["controlCoeff"]
- config["bPeaks"]["log2FC"]
- config["bPeaks"]["averageQuantiles"]
- config["bPeaks"]["windowSize"]
- config["bPeaks"]["windowOverlap"]

bwa_index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rule for the creation of BWA index. Has to be done only once. 
The output file is used to test whether the index already exists when aligning.

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["version"]
- config["genome"]["fasta_file"]

bwa
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Requires the indexing to have previously been done (using the rule bwa_index).

Required parameters:

- config["dir"]["fastq"]
- config["dir"]["samples"]
- config["genome"]["version"]
- config["genome"]["fasta_file"]
- config["qsub"]

Optional parameters:

- config["bwa"]["dir"]
- config["bwa"]["threads"]

count_reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A set of rules to count the number of reads in NGS files with different formats.

Includes: 

- rule count_reads_fastq

Count number of reads in a fastq-formatted file (unaligned reads).

- rule count_reads_fastq_gz

Count number of reads in a gzipped fastq-formatted file (unaligned reads).

- rule count_reads_bam

Count number of reads in a bam-formatted file (binary alignment map, compressed sam).

 -rule count_reads_sam

Count number of reads in a bam-formatted file (binary alignment map, compressed sam).

 -rule count_features_bed

Count number of features in a bed-formatted file.


cufflinks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Annotate a genome (infer transcript regions) based on RNA-seq data.
Assemble transcripts from RNA-seq data and produce a file with the location of detected transcripts.

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["gtf_file"]

Optional parameters:

- config["cufflinks"]["libtype"]
- config["cufflinks"]["threads"]

cutadapt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Trimming of raw read files using the tool cutadapt and the wrapper Trim Galore. 

Currently only working with single end data. 

Required parameters:

- config["qsub"]
- config["metadata"]["seq_type"]

Optional parameters:

- config["cutadapt"]["length_threshold"]
- config["cutadapt"]["qual_threshold"]
- config["metadata"]["strands"]


dot_graph
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This rule generates dot files for snakemake's DAG and rulegraph. 

Required parameter: 

- config["metadata"]["configfile"]

dot_to_image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Following rule dot_graph, this rule creates png, pdf and svg outputs for dot graphs from snakemake. 

fastqc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Check the quality of the reads in a fastq or a bam file using the program fastQC (quality control). 
Results are stored in folder named '{reads}_fastqc'.

Custom parameters specified in the configuration file with the option config["fastqc"]["other options"] will be passed to fastQC. 

Required parameters:

- config['qsub']

Optional parameters:

- config['fastqc']['other_options']

featnb_from_bed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Count number of features in a bed-formatted file.

Required parameters:

- config["qsub"]

genome_coverage_bedgraph_strands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute two strand-specific genome coverage files from a bam-formatted file with aligned reads. 
The coverage files are in bedgraph format, which can be loaded in the genome viewer IGV.

Required parameters: 

- config["qsub"]


genome_coverage_bedgraph
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute genome coverage from a bam-formatted file with aligned reads. 
The coverage file is in bedgraph format (with extension .bedgraph), which can be loaded in the genome viewer IGV. 

Note however that IGV issues a warning when bedgraph files are given in input, and recommends to use the tdf format instead. 
We implemented hereafter rules to convert bedgraph to tdf.

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["fasta_file"]

genome_coverage_bigwig_normalized
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Uses bamCompare tool from the deepTools suite.

Required parameters:

- config["qsub"]

genome_coverage_bigwig
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute genome coverage from a bam-formatted file with aligned reads and produce a bigWig file. 
Uses bamCoverage tool from the deepTools suite. 

Required parameters:

- config["qsub"]


genome_coverage_dz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute coverage (reads per position) for each position of a genome, from a bam-formatted file with aligned reads.

BEWARE: this rule is useful for small genomes (Bacteria, Fungi) but would produce a very big file for Metazoan or Plant genomes.

Required parameters:

- config["qsub"]

genome_coverage_wig
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute genome coverage from a bam-formatted file with aligned reads and produce a wig file, the recommended format to upload coverage-type data as UCSC tracks.

Required parameters:

- config["qsub"]

get_chrom_sizes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This rule generates a file containg the chromosome sizes, with file extension *.genome. 
This file is required by a number of bedtools utilities.

Required parameters:

- config["qsub"]

getfasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Get fasta from bed file using the bedtools suite.

Required parameters:

- config["qsub"]
- config["dir"]["genome"] 
- config["genome"]["fasta_file"]

gunzip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Uncompress a file with the gunzip program. 
The rule is very simple, but is convenient to use in a workflow: 
it can be used to fix some dependencies on.gz extensions, and/or to send compression jobs to a queue.

Required parameters:

- config["qsub"]

gzip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Uncompress a file with the gunzip program. 
The rule is very simple, but is convenient to use in a workflow: 
it can be used to fix some dependencies on.gz extensions, and/or to send compression jobs to a queue.

Required parameters:

- config["qsub"]


homer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peak-calling with HOMER software, findPeaks algorithm.
Input formats: .sam, .bam, .bed (bam input requires samtools to be installed)

The genome parameter can be either 
the code of a genome installed in Homer (eg HG18, dm3...) or 
a fasta file (see http://homer.salk.edu/homer-fdr{fdr}_peaks/introduction/update.html)

Required parameters:

- config["qsub"]
- config["dir"]["samples"]
- config["dir"]["peaks"]
- config["genome"]["fasta_file"]
- config["genome"]["size"]

Optional parameters:

- config["homer"]["style"]
- config["homer"]["L"]
- config["homer"]["F"]
- config["homer"]["P"]
- config["homer"]["fdr"]

index_bam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Index a bam file by creating a .bai file with Samtools
The input bam MUST be sorted by position (rule bam_by_pos). 

Required parameters:

- config["qsub"]

index_fasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Index a fasta file by creating an .fai file with Samtools.

Required parameters:

- config["qsub"]


macs14
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peak-calling with macs14.
Input: bed (others supported)
Output: bed

Required parameters:

- config["genome"]["size"]
- config["qsub"]
- config["dir"]["samples"]
- config["dir"]["peaks"]

Optional parameters:

- config["macs14"]["pval"]
- config["macs14"]["keep_dup"]
- config["macs14"]["bandwidth"]
- config["macs14"]["mfold"]
- config["macs14"]["other_options"]


macs2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peak-calling with MACS2.
Input: bed
Output: bed

Required parameters:

- config["dir"]["samples"]
- config["dir"]["peaks"]
- config["genome"]["size"]
- config["qsub"]

Optional parameters:

- config["macs2"]["qval"]
- config["macs2"]["keep_dup"]
- config["macs2"]["band_width"]
- config["macs2"]["mfold_min"]
- config["macs2"]["mfold_max"]
- config["macs2"]["other_options"]
- config["macs2"]["type"]


matrix_clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*This rule is currently incomplete.*

Motif discovery using the peak-motifs pipeline.

Documentation of tool `here <http://floresta.eead.csic.es/rsat/help.peak-motifs.html>`__.

Required parameters:

- config["qsub"]

Optional parameters:

matrix_quality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Measuring peak enrichment for motifs. 
Requires at least 2 sets of peaks and 1 motif database. 

Documentation of tool `here <http://embnet.ccg.unam.mx/rsa-tools/help.matrix-quality.html>`__.

Required parameters:

- config["qsub"]
- config["matrix-quality"]["background"]

Optional parameters:


md5sum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compute the md5sum signature for a given file, which enables to check the consistency of its content after transfer. 

Note: md5sum is recommended for submitting NGS data to GEO.

Usage: integrate in the targets of a workflow. 
Alternatively, can be called directly on the command line with find.

Example: find all the fastq files in a directory (named fastq) and 
compute one md5sum file for each, and assign the takss to 20 jobs in the scheduler.

::

    find fastq/ -name '*.fastq'  \
       | awk '{print $1".md5sum"}' \
       | xargs snakemake  -j 20 -p \
           -s gene-regulation/scripts/snakefiles/rules/md5sum.rules \
           --configfile metadata/Glossina_palpalis.yml

Required parameters:

- config["qsub"]

merge_lanes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Concatenate multiple fastq files to produce a merged fastq file.

This rule typically serves to merge raw reads (fastq) corresponding to
multiple sequencing lanes for the same sample into a single fastq file
per sample.

Since the file naming conventions are highly dependent on the sequencing
platform, the file grouping is read from a user-provided text file with
tab-separated values (extension .tsv). This file must have been specified
in the config file, as config["metadata"]["lane_merging"].

This file must contain at least two columns with this precise header:

- source_file
- merged_file

Additional columns can be included but will be ignored.

There should be a N to 1 correspondence from source file to merge file
  (each source file should in principle be assigned to a single merged file).

Source files are supposed to be compressed fastq sequence files (.fastq.gz).

The output file is an uncompressed fastq file, because bowtie version 1 does not support gzipped files as input.

Required parameters:

- config["metadata"]["lane_merging"]  file indicating the source/merged file names
- config["dir"]["fastq"]              base of the directory containing the fastq files


mosaics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peak-calling with mosaics R package.
Input: bed
Output: bed

Required parameters:

- config["dir"]["samples"]
- config["qsub"]

Optional parameters:

- config["mosaics"]["FDR"]
- config["mosaics"]["frag_len"]
- config["mosaics"]["bin_size"]
- config["mosaics"]["type"]

peak_motifs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Motif discovery using the peak-motifs pipeline from `RSAT <rsat.eu>`__.

`Documentation <http://floresta.eead.csic.es/rsat/help.peak-motifs.html>`__.

Required parameters:

- config["qsub"]
- config["peak-motifs"]["motif_db"]

Optional parameters:

- config["peak-motifs"]["tasks"]
- config["peak-motifs"]["disco"]


readnb_from_bam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Count number of reads in a bam-formatted file (binary alignment map, compressed sam).

Required parameters:

- config["qsub"]

readnb_from_fastq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Count number of reads in a fastq-formatted file (unaligned reads).

Required parameters:

- config["qsub"]

readnb_from_sam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Count number of reads in a bam-formatted file (binary alignment map, compressed sam).

Required parameters:

- config["qsub"]

sam_to_bam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Convert reads from SAM (sequence alignment map) to BAM (binary alignment map) format.

Required parameters:

- config["qsub"]

sartools_DESeq2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This rule is designed to perform differential expression analysis of RNA-seq data with DESeq2, using the R package `SARTools <https://github.com/PF2-pasteur-fr/SARTools/>`__. 

It requires replicated data to run properly.

Required parameters:

- config["qsub"]
- config["author"]
- config["dir"]["samples"]
- config["dir"]["diffexpr"]

Optional parameters:

- config["DESeq2"]["featuresToRemove"]
- config["DESeq2"]["varInt"]
- config["DESeq2"]["condRef"]
- config["DESeq2"]["batch"]
- config["DESeq2"]["alpha"]
- config["DESeq2"]["pAdjustMethod"]
- config["DESeq2"]["fitType"]
- config["DESeq2"]["cooksCutoff"]
- config["DESeq2"]["independentFiltering"]
- config["DESeq2"]["typeTrans"]
- config["DESeq2"]["locfunc"]


sartools_edgeR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This rule is designed to perform differential expression analysis of RNA-seq data with edgeR, using the R package `SARTools <https://github.com/PF2-pasteur-fr/SARTools>`__.

It requires replicated data to run properly.

Required parameters:

- config["qsub"]
- config["author"]
- config["dir"]["samples"]
- config["dir"]["diffexpr"]

Optional parameters:

- config["edgeR"]["featuresToRemove"]
- config["edgeR"]["varInt"]
- config["edgeR"]["condRef"]
- config["edgeR"]["batch"]
- config["edgeR"]["alpha"]
- config["edgeR"]["pAdjustMethod"]
- config["edgeR"]["fitType"]
- config["edgeR"]["cpmCutoff"]
- config["edgeR"]["gene_selection"]
- config["edgeR"]["normalizationMethod"]

sartools_targetfile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This rule creates a so-called "targetfile", which is required by SARTools to run differential expression analyses with rules sartools_DESeq2 and sartools_edgeR.

Required parameters:

- config["qsub"]
- config["dir"]["samples"]
- config["dir"]["diffexpr"]

sickle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Trimming raw reads with sickle.

Required parameters:

- config["qsub"]
- config["metadata"]["seq_type"]

Optional parameters:

- config["sickle"]["qual_threshold"]
- config["sickle"]["length_threshold"]
- config["sickle"]["qual_type"]
- config["metadata"]["strands"]

split_bam_by_strands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Split a bam file in two separate bam files containing respectively the reads on the plus and minus strand. 

Required parameters:

- config["qsub"]

spp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peak-calling with SPP (R package).
Input: bam
Output: narrowPeak + bed format

Required parameters:

- config["qsub"]
- config["dir"]["samples"]
- config["dir"]["peaks"]

Optional parameters:

- config["spp"]["fdr"]


sra_to_fastq_split
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Converts SRA files in separate fastq files for paired-end data with SRA toolkit. 

Required parameters:

- config["qsub"]
- config["dir"]["reads_source"]
- config["dir"]["fastq"]

Optional parameters:

- config["fastq_dump"]["options"]

Usage example :

IMPORT = expand(FASTQ_DIR + "/{samples}/{samples}.fastq", samples=SAMPLE_IDS) 

sra_to_fastq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Converts SRA files in fastq format with SRA toolkit (songle-end data only). 

Required parameters:

- config["qsub"]
- config["dir"]["reads_source"]
- config["dir"]["fastq"]


subread_align
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Align each sample with the R-package subread.

To align each sample on the reference genome the R-package subread first needs to build a index with the function builindex(). 
The alignment is then executed with the function align(), which calls the tool read mapping tool Subread.  

Reference: Liao Y, Smyth GK and Shi W (2013). The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. 
Nucleic Acids Research, 41(10):e108

Required parameters:

- config["dir"]["fastq"]
- config["dir"]["samples"]
- config["genome"]["version"]
- config["genome"]["fasta_file"]
- config["qsub"]
- config["metadata"]["seq_type"]

Optional parameters:

- config["subread-align"]["dir"]
- config["subread-align"]["threads"]
- config["subread-align"]["max_mismatches"]
- config["subread-align"]["align_options"]
- config["subread-align"]["seq_data"]


subread_featureCounts_all
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This rule computes bam files and produces a tab file containing feature counts for all samples, 
using featureCounts from the subread toolkit. 

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["gtf_file"]

Optional parameters:

- config["subread-featureCounts"]["attr_type"]
- config["subread-featureCounts"]["feature_type"]
- config["subread-featureCounts"]["multi_mapping"]
- config["subread-featureCounts"]["strand_specificity"]

Usage: 

::

    featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 


subread_featureCounts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


This rule computes bam files and produces tab files containing feature counts for every sample separately, 
using featureCounts from the subread toolkit. 

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["gtf_file"]

Optional parameters:

- config["subread-featureCounts"]["attr_type"]
- config["subread-featureCounts"]["feature_type"]
- config["subread-featureCounts"]["multi_mapping"]
- config["subread-featureCounts"]["strand_specificity"]

Usage: 

::

    featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 

subread_index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rule for the creation of subread index. Has to be done only once.
The output file is used to test whether the index already exists when aligning.

Reference: Liao Y, Smyth GK and Shi W (2013). The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. 
Nucleic Acids Research, 41(10):e108

Required parameters:

- config["qsub"]
- config["dir"]["genome"]
- config["genome"]["version"]
- config["genome"]["fasta_file"]



swembl
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Peak-calling with SWEMBL.

Beware: for SWEMBL the peaks MUST be sorted by position, other wise SWEMBL runs indefinitely. 
Usually by default we sort all bam files by position after alignment (with rule bam_by_pos). 

Required parameters:

- config["qsub"]
- config["dir"]["samples"]
- config["dir"]["peaks"]

Optional parameters:

- config["swembl"]["x"]
- config["swembl"]["R"]
- config["swembl"]["N"]

tophat
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read mapping for single or paired end data using Tophat. 
Requires the indexing to have previously been done (using the rule bowtie2_index).

Required parameters:

- config["dir"]["fastq"]
- config["dir"]["samples"]
- config["metadata"]["seq_type"]
- config["dir"]["genome"]
- config["genome"]["version"]
- config["genome"]["fasta_file"]
- config["qsub"]

Optional parameters:

- config["tophat"]["max_mismatches"]
- config["tophat"]["threads"]



Python scripts (.py)
----------------------------------------------------------------

*todo*

R scripts
----------------------------------------------------------------

*todo*


Configuration files (yaml)
----------------------------------------------------------------

*todo*


R markdown files (.Rmd)
----------------------------------------------------------------

*todo*

Tabulated files (.tab)
----------------------------------------------------------------

We use tabulated files in order to define and describe the samples 
to be processed in the workflows. 

Examples of these files are available in the *examples* folder of the 
library. 

Sample description files (samples.tab)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*todo*

Experimental design files (design.tab)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*todo*
