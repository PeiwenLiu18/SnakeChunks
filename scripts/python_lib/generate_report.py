## Draft python script to generate reports

# Usage: python generate_report.py path_to_configfile

import os
import sys
import pandas as pd
import yaml
import collections

import web
import util
import statistics

dataset_id = snakemake.params["dataset"]
outputdir = snakemake.params["outdir"]
project_name = snakemake.params["project_name"]

# Framework/tools dictionary

#...

# Samples dictionary (separate ChIP and RNA?)

samples_ids = snakemake.params["samples_ids"]
samples_dir = snakemake.params["samples_dir"]
trim_tools = snakemake.params["trim_tools"]

d_samples = {}

for sample in samples_ids:
    path = samples_dir + "/" + sample + "/" + sample + ".fastq"
    d_samples[sample] = {}
    d_samples[sample]["path"] = path



if __name__ == '__main__':
    web.indexPage(outputdir, project_name, dataset_id)
    web.homePage(outputdir, project_name, dataset_id, "home", description="")
    web.writeCategory(outputdir, project_name, dataset_id, "chip", "chip", d_samples)
    web.downloadPage(outputdir, project_name, dataset_id, "download", description="")




### Read config.yml
##   analysis steps
##   filenames

#if not ("tools" in config.keys()):
#    sys.exit("The parameter config['tools'] should be specified in the config file.")

#step_list = []

## Check if there's a trimming step
#if ("trimming" in config["tools"].keys()):
#    trimming_step = True
#    step_list.append("trimming")
#    trimming_tools = config["tools"]["trimming"].split()

#else:
#    trimming_step=False

## Check if there's a mapping step
#if ("mapping" in config["tools"].keys()):
#    mapping_step = True
#    step_list.append("trimming")
#    mapping_tools = config["tools"]["mapping"].split()

#else:
#    mapping_step=False


## Check if there's a peak-calling step
#if ("peakcalling" in config["tools"].keys()):
#    peakcalling_step = True
#    step_list.append("trimming")
#    peakcalling_tools = config["tools"]["peakcalling"].split()

#else:
#    peakcalling_step = False




#if TRIMMING_TOOLS:
#    PREFIX = expand("{trimmer}_{aligner}", aligner=MAPPING_TOOLS, trimmer=TRIMMING_TOOLS)
#else:
#    PREFIX = expand("{aligner}", aligner=MAPPING_TOOLS)


### Read files -> script python

#    # If fastq files
#        # Generate Fastq stats
#        # Write to dictionary
#        # Copy files to reports dir

#    # If bam files
#        # Generate bam stats
#        # Write to dictionary
#        # Copy tdf files to reports dir

#    # If bed files
#        # Generate bed stats
#        # Write to dictionary
#        # Copy files to reports dir


## Write report



    # Write header

        # If fastq files
            # Write block

        # If bam files
            # Write block

        # If bed files
            # Write block

        # Plots (todo)

    # Write footer

## Write genome viewer session

    # IGV

    # Genome viewer
