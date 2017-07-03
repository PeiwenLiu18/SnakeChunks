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

Sort aligned reads (in bam format) by positions using 'samtools sort'. 

Required parameters:

- config["qsub"]

bam_by_pos
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sort aligned reads (in bam format) by name, using 'samtools sort'. 

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

bwa
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

count_reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cufflinks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cutadapt
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

dot_graph
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

dot_to_image
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

fastqc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

featnb_from_bed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

genome_coverage_bedgraph_strands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

genome_coverage_bedgraph
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

genome_coverage_bigwig_normalized
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

genome_coverage_bigwig
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

genome_coverage_dz
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

genome_coverage_wig
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get_chrom_sizes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

getfasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

gunzip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

gzip
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

homer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

index_bam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

index_fasta
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

macs14
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

macs2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

matrix_clustering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

matrix_quality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

md5sum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

merge_lanes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

mosaics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

readnb_from_fastq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

readnb_from_sam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sam_to_bam
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sartools_DESeq2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sartools_edgeR
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sartools_targetfile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sickle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

split_bam_by_strands
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

spp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sra_to_fastq_split
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sra_to_fastq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

subread_align
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

subread_featureCounts_all
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

subread_featureCounts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

subread_index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

swembl
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

tophat
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




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
