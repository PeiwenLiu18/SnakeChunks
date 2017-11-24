################################################################
## This makefile downloads data for tutorials based on dataset GSE41195
## 
## Ref:     Myers KS, Yan H, Ong IM, Chung D et al. Genome-scale analysis of escherichia coli FNR reveals complex features of transcription factor binding. 
## 			PLoS Genet 2013 Jun;9(6):e1003565. PMID: 23818864
##
## Author: Claire Rioualen
## Date: 2017-11-24
##
## Usage: make -f tutorial_material.mk all ANALYSIS_DIR=/my/dir
##
export LC_ALL=C
export LANG=C

.PHONY: \
	init\
	create_dir\
	download_genome_data\
	download_raw_data



create_dir:
	@echo "Creating ANALYSIS_DIR directory" && \
	export ANALYSIS_DIR=$(ANALYSIS_DIR) && \
	mkdir -p $(ANALYSIS_DIR) && \
	ln -sf $(HOME)/SnakeChunks-4.0 $(ANALYSIS_DIR)/SnakeChunks


### Download genome & annotations 


download_genome_data:
	mkdir -p ${ANALYSIS_DIR}/genome && \
	wget -nc ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz -P ${ANALYSIS_DIR}/genome && \
	wget -nc ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.chromosome.Chromosome.gff3.gz -P ${ANALYSIS_DIR}/genome && \
	gunzip ${ANALYSIS_DIR}/genome/*.gz


### Download ChIP-seq data 

download_raw_data:
	cd $(ANALYSIS_DIR) && \
	mkdir -p ${ANALYSIS_DIR}/ChIP-seq/fastq/GSM1010224 ${ANALYSIS_DIR}/ChIP-seq/fastq/GSM1010220 && \
	wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576934/SRR576934.fastq.gz -P ${ANALYSIS_DIR}/ChIP-seq/fastq/GSM1010220 && \
	wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR576/SRR576938/SRR576938.fastq.gz -P ${ANALYSIS_DIR}/ChIP-seq/fastq/GSM1010224 && \
	gunzip -c ${ANALYSIS_DIR}/fastq/GSM1010220/SRR576934.fastq.gz > ${ANALYSIS_DIR}/fastq/GSM1010220/GSM1010220.fastq; rm -f ${ANALYSIS_DIR}/fastq/GSM1010220/SRR576934.fastq.gz && \
	gunzip -c ${ANALYSIS_DIR}/fastq/GSM1010224/SRR576938.fastq.gz > ${ANALYSIS_DIR}/fastq/GSM1010224/GSM1010224.fastq; rm -f ${ANALYSIS_DIR}/fastq/GSM1010224/SRR576938.fastq.gz


### Copy metadata from SnakeChunks library to analysis directory

copy_metadata:
	mkdir metadata ; cp SnakeChunks/examples/GSE41187/* metadata

### all

all: \
	init\
	create_dir\
	copy_metadata\
	download_genome_data\
	download_raw_data

