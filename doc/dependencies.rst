Dependencies
================================================================

*Note: this section needs to be refreshed*

Manual installation
----------------------------------------------------------------

This manual aims at helping you install the necessary programs and
dependencies in order to have the snakemake workflows work. It was
designed for Unix-running computers (Ubuntu, Debian).

General requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generic tools
****************************************************************

ssh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    sudo apt-get install ssh
    ssh-keygen

rsync
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`rsync <https://rsync.samba.org/>`__ is an open source utility that
provides fast incremental file transfer.

::

    sudo apt-get install rsync

git
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Install git on your machine.

::

    sudo apt-get install git

Optional:

-  Create an account on `GitHub <https://github.com>`__.
-  Add your ssh public key to your GitHub account settings (account >
   settings > SSH keys > add SSH key).

::

    less ~/.ssh/id_rsa.pub

zlib
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several tools require this dependency (e.g. sickle, bamtools...).

::

    sudo apt-get install libz-dev

qsub
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



Create bin/ and app\_sources/ (optional)
****************************************************************

While some programs will be installed completely automatically, others 
will not. Here we create a directory that will be used for manual
installations.

::

    mkdir $HOME/bin
    mkdir $HOME/app_sources

You might then have to edit your ``$PATH`` manually (see next section).

Edit ``$PATH``
****************************************************************

In order to use manually installed programs and make them executable,
you may have to update your ``$PATH`` environment variable. You can do
so by editing the ``~/.profile`` file.

::

    nano ~/.profile

Fetch this paragraph and add the path to manually installed executables:

::

    # set PATH so it includes user's private bin if it exists
    if [ -d "$HOME/bin" ] ; then
        PATH="$HOME/bin:$PATH"
    fi

Execute the file to validate the change.

::

    source ~/.profile

Snakemake workflows basic requirements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python
****************************************************************

Snakemake requires to have Python 3.3+ installed. 
You can check this by issuing the following commands in a terminal:

::

    python --version # usually the default python version is 2.7+
    python3 --version

If you don't have python 3 you should install it.

::

    sudo apt-get install python3

Install pip and pip3.

::

    sudo apt-get install python-pip
    sudo apt-get install python3-pip

Not installed natively?

::

    apt-get install python-dev
    apt-get install python3.4-dev

Pandas library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This library is used in order to read tab-delimited files used in the workflows 
(see files ``samples.tab`` and ``design.tab``).

::

    pip3 install pandas

Package rpy2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    pip3 install "rpy2<2.3.10"



R
****************************************************************


*todo*

Snakemake
****************************************************************

-  `Documentation <http://snakemake.readthedocs.io>`__
-  `FAQ <https://bitbucket.org/snakemake/snakemake/wiki/FAQ>`__
-  `Forum <https://groups.google.com/forum/#!forum/snakemake>`__
-  See also Snakemake section for tutorials. 

Now you have installed Python 3 and pip3 (see previous section), you can
install snakemake safely.

::

    pip3 install snakemake

You can check that snakemake works properly with this basic script:

::

    """Snakefile to test basic functions of snakemake.
    """
    rule all:
        input: expand("bye.txt")

    rule hello:
        """Write HELLO in a text file named hello.txt.
        """
        output: "hello.txt"
        message: "Generating {output} file."
        shell: "echo HELLO > {output}"

    rule bye:
        """Write BYE in a text file named bye.txt.
        """
        input: "hello.txt"
        output: "bye.txt"
        message: "Generating {output} file."
        shell: "echo BYE > {output}"

-  Save it to ``~/workspace/hello.py``.
-  Issue the command ``cd ~/workspace ; snakemake -s hello.py``.
-  2 files should be created: ``hello.txt`` and ``bye.txt``.

As of December 2015, you need snakemake version 3.4+.

::

    pip3 install snakemake --upgrade

If you want to use Snakemake reports function (optional):

::

    pip3 install docutils

Graphviz
****************************************************************

Snakemake can generate useful graphviz outputs.

::

    sudo apt-get install graphviz

NGS analysis software & tools
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Quality assessment
****************************************************************

FastQC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`FastQC <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>`__
aims to provide a simple way to do some quality control checks on raw
sequence data coming from high throughput sequencing pipelines. It
provides a modular set of analyses which you can use to give a quick
impression of whether your data has any problems of which you should be
aware before doing any further analysis.

Links:

-  `QC Fail Sequencing <https://sequencing.qcfail.com/>`__

-  `FastQC results
   interpretation <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/>`__

FastQC is available from linux repositories:

::

    sudo apt-get install fastqc

However, since it's an older version, it can cause problems of dependencies. 

We recommend installing it manually: 

::

    cd $HOME/app_sources
    wget --no-clobber http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
    unzip -o fastqc_v0.11.5.zip
    chmod +x FastQC/fastqc
    ln -s -f $HOME/app_sources/FastQC/fastqc $HOME/bin/fastqc

Trimming
****************************************************************

Sickle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Sickle <https://github.com/najoshi/sickle>`__ is a trimming tool which
better the quality of NGS reads.

-  Pre-requisite: install ``zlib`` (*link to section*).
-  Clone the git repository into your bin (*link to section*) and run
   ``make``.

::

    cd $HOME/app_sources
    git clone https://github.com/najoshi/sickle.git 
    cd sickle 
    make 
    cp sickle $HOME/bin



Alignment/mapping
****************************************************************

BWA
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`BWA <http://bio-bwa.sourceforge.net/>`__ is a software package for
mapping low-divergent sequences against a large reference genome, such
as the human genome.

-  `Manual <http://bio-bwa.sourceforge.net/bwa.shtml>`__

-  `Publication <http://www.ncbi.nlm.nih.gov/pubmed/19451168>`__ 

Li H. and Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.

::

    sudo apt-get install bwa

Bowtie
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    cd $HOME/app_sources
    wget --no-clobber http://downloads.sourceforge.net/project/bowtie-bio/bowtie/$(BOWTIE1_VER)/bowtie-$(BOWTIE1_VER)-linux-x86_64.zip
    unzip bowtie-$(BOWTIE1_VER)-linux-x86_64.zip
    cp `find bowtie-$(BOWTIE1_VER)/ -maxdepth 1 -executable -type f` $HOME/bin


Bowtie2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`General
documentation <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`__

`Instructions <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2>`__

`Downloads <https://sourceforge.net/projects/bowtie-bio/files/bowtie2/>`__

::

    cd $HOME/app_sources
    wget http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip
    unzip bowtie2-2.2.6-linux-x86_64.zip
    p `find bowtie2-$(BOWTIE2_VER)/ -maxdepth 1 -executable -type f` $HOME/bin


Peak-calling
****************************************************************

bPeaks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Web page <http://bpeaks.gene-networks.net/>`__

Peak-caller developped specifically for yeast, can be useful in order to
process small genomes only.

Available as an R library.

::

    install.packages("bPeaks")
    library(bPeaks)

HOMER
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Web page <http://homer.salk.edu/>`__

`Install
instructions <http://homer.salk.edu/homer/introduction/install.html>`__

::

    mkdir $HOME/app_sources/homer
    cd $HOME/app_sources/homer
    wget "http://homer.salk.edu/homer/configureHomer.pl"
    perl configureHomer.pl -install homer
    cp `find $HOME/app_sources/homer/bin -maxdepth 1 -executable -type f` $HOME/bin

The basic Homer installation does not contain any sequence data. To
download sequences for use with HOMER, use the configureHomer.pl script.
To get a list of available packages:

::

    perl $HOME/bin/HOMER/configureHomer.pl -list

To install packages, simply use the -install option and the name(s) of
the package(s).

::

    perl  $HOME/bin/HOMER/configureHomer.pl -install mouse # (to download the mouse promoter set)
    perl  $HOME/bin/HOMER/configureHomer.pl -install mm8   # (to download the mm8 version of the mouse genome)
    perl  $HOME/bin/HOMER/configureHomer.pl -install hg19  # (to download the hg19 version of the human genome)

Supported organisms:

+-----------------+--------------------+
| Organism        | Assembly           |
+=================+====================+
| Human           | hg17, hg18, hg19   |
+-----------------+--------------------+
| Mouse           | mm8, mm9, mm10     |
+-----------------+--------------------+
| Rat             | rn4, rn5           |
+-----------------+--------------------+
| Frog            | xenTro2, xenTro3   |
+-----------------+--------------------+
| Zebrafish       | danRer7            |
+-----------------+--------------------+
| Drosophila      | dm3                |
+-----------------+--------------------+
| C. elegans      | ce6, ce10          |
+-----------------+--------------------+
| S. cerevisiae   | sacCer2, sacCer3   |
+-----------------+--------------------+
| S. pombe        | ASM294v1           |
+-----------------+--------------------+
| Arabidopsis     | tair10             |
+-----------------+--------------------+
| Rice            | msu6               |
+-----------------+--------------------+

HOMER can also work with custom genomes in FASTA format and gene
annotations in GTF format.

MACS 1.4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  `Documentation <http://liulab.dfci.harvard.edu/MACS/00README.html>`__
-  `Installation manual <http://liulab.dfci.harvard.edu/MACS/INSTALL.html>`__

::

    cd $HOME/app_sources
    wget "https://github.com/downloads/taoliu/MACS/MACS-1.4.3.tar.gz"
    tar -xvzf MACS-1.4.3.tar.gz
    cd MACS-1.4.3
    sudo python setup.py install
    macs14 --version


MACS2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  `Webpage <https://github.com/taoliu/MACS/>`__

::

    sudo apt-get install python-numpy
    sudo pip install MACS2


SPP R package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This one might be a little but tricky (euphemism).

Several possibilities, none of which have I had the courage to retry lately. 

- In R

::

    source("http://bioconductor.org/biocLite.R")
    biocLite("spp")
    install.packages("caTools")
    install.packages("spp")

- In commandline

::

    apt-get install libboost-all-dev
    cd $HOME/app_sources
    wget -nc http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_1.11.tar.gz
    sudo R CMD INSTALL spp_1.11.tar.gz

- Using git (I haven't tried this one but it looks more recent) (see `github page <https://github.com/hms-dbmi/spp>`__)

::

    require(devtools)
    devtools::install_github('hms-dbmi/spp', build_vignettes = FALSE)


I also wrote a little protocol a while ago. 
Here's the procedure on Ubuntu 14.04, in this very order:

In unix shell:

::

    # unix libraries
    apt-get update
    apt-get -y install r-base
    apt-get -y install libboost-dev zlibc zlib1g-dev

In R shell:

::

    # Rsamtools
    source("http://bioconductor.org/biocLite.R")
    biocLite("Rsamtools")

In unix shell:

::

    # spp
    wget http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_1.11.tar.gz
    sudo R CMD INSTALL spp_1.11.tar.gz

A few links:

-  Download page can be found
   `here <http://compbio.med.harvard.edu/Supplements/ChIP-seq/>`__,
   better chose version ``1.11``.
-  SPP requires the Bioconductor library
   `Rsamtools <https://bioconductor.org/packages/release/bioc/html/Rsamtools.html>`__
   to be installed beforehand.
-  Unix packages ``gcc`` and ``libboost`` (or equivalents) must be
   installed.
-  You can find a few more notes
   `here <http://seqanswers.com/forums/archive/index.php/t-22653.html>`__.
-  Good luck!

SWEMBL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  `SWEMBL beginner's
   manual <http://www.ebi.ac.uk/~swilder/SWEMBL/beginners.html>`__

::

    cd $HOME/app_sources
    wget "http://www.ebi.ac.uk/~swilder/SWEMBL/SWEMBL.3.3.1.tar.bz2"
    bunzip2 -f SWEMBL.3.3.1.tar.bz2
    tar xvf SWEMBL.3.3.1.tar
    rm SWEMBL.3.3.1.tar
    chown -R ubuntu-user SWEMBL.3.3.1
    cd SWEMBL.3.3.1
    make

It seems there could be issues with C flags. To be investigated. 

Motif discovery, motif analysis
****************************************************************

RSAT suite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*see dedicated section*

Miscellaneous
****************************************************************

SRA toolkit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This toolkit includes a number of programs, allowing the conversion of
``*.sra`` files. ``fastq-dump`` translates ``*.sra`` files to
``*.fastq`` files.

-  `SRA format <http://www.ncbi.nlm.nih.gov/Traces/sra/>`__
-  `fastq-dump
   manual <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump>`__
-  `Installation
   manual <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std>`__

You can download last version
`here <http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>`__,
or issue the following commands:

::

    cd $HOME/app_sources
    wget -nc http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz
    tar xzf sratoolkit.2.5.2-ubuntu64.tar.gz
    cp `find sratoolkit.2.5.2-ubuntu64/bin -maxdepth 1 -executable -type l` $HOME/bin

You can also install SRA toolkit simply by issuing this
command, but likely it won't be the most recent release:

::

    sudo apt-get install sra-toolkit

::

    fastq-dump --version
      fastq-dump : 2.1.7

Samtools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SAM (Sequence Alignment/Map) format is a generic format for storing
large nucleotide sequence alignments.

`SAMtools <http://samtools.sourceforge.net/>`__ provides several tools
to process such files.

::

    cd $HOME/app_sources
    wget --no-clobber http://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2
    bunzip2 -f samtools-1.3.tar.bz2
    tar xvf samtools-1.3.tar
    cd samtools-1.3
    make 
    sudo make install

Bedtools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `bedtools <http://bedtools.readthedocs.org/en/latest/>`__ utilities
are a swiss-army knife of tools for a wide-range of genomics analysis
tasks. For example, bedtools allows one to intersect, merge, count,
complement, and shuffle genomic intervals from multiple files in
widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF.

::

    sudo apt-get install bedtools

or get the latest version:

::

    cd $HOME/app_sources
    wget --no-clobber https://github.com/arq5x/bedtools2/releases/download/v2.24.0/bedtools-2.24.0.tar.gz
    tar xvfz bedtools-2.24.0.tar.gz
    cd bedtools2
    make
    sudo make install



Bedops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    cd $HOME/app_sources
    wget -nc https://github.com/bedops/bedops/releases/download/v2.4.19/bedops_linux_x86_64-v2.4.19.tar.bz2
    tar jxvf bedops_linux_x86_64-v2.4.19.tar.bz2
    mkdir bedops
    mv bin bedops
    cp bedops/bin/* $HOME/bin

Deeptools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    cd $HOME/app_sources
    git clone https://github.com/fidelram/deepTools
    cd deepTools
    python setup.py install

Picard tools 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*todo*



Makefile
----------------------------------------------------------------

*Has to be revised*

The Gene-regulation library comprises a makefile that can install most of the 
dependencies described in the previous section. 

It currently allows running the following workflows:

- import_from_sra.wf
- quality_control.wf
- ChIP-seq.wf

::

    cd $GENE_REG_PATH
    make -f gene-regulation/scripts/makefiles/install_tools_and_libs.mk all
    source ~/.bashrc

Conda
----------------------------------------------------------------

*This section has to be written*

::

    conda install -c bioconda sickle=0.5 
    conda install -c bioconda bowtie=1.2.0 
    conda install -c bioconda bowtie2=2.3.0 
    conda install -c bioconda subread=1.5.0.post3 
    conda install -c bioconda tophat=2.1.1 
    conda install -c bioconda bwa=0.7.15 
    conda install -c bioconda fastqc=0.11.5 
    conda install -c bioconda macs2=2.1.1.20160309 
    conda install -c bioconda homer=4.8.3 
    conda install -c bioconda bedtools=2.26.0 
    conda install -c bioconda samtools=1.3.1 
    conda install -c bioconda bamtools=2.4.0 


