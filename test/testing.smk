#!/usr/bin/env snakemake -s
import os
from snakemake.utils import min_version
# set path for configuration file
configfile: "/home/glbrc.org/millican/repos/trait-mapper/test/config.yml"
# set working directory path
workdir: "/home/glbrc.org/millican/repos/trait-mapper/workflow"

def get_sample_names():
    sample_list = []
    for sample in config["input_samples"]:
        bname = os.path.basename(sample)
        name = bname.replace(".fastq.gz", "")
        sample_list.append(name)

# samples = get_sample_names() # this is not for the testing phase

module search_workflow:
    snakefile:
        "/home/glbrc.org/millican/repos/trait-mapper/workflow/search_workflow.smk"
    config:
        config

use rule * from search_workflow as search_*

# import modules
module taxonomy_workflow:
    snakefile:
        "/home/glbrc.org/millican/repos/trait-mapper/workflow/taxonomy_workflow.smk"
    config:
        config
        
use rule * from taxonomy_workflow as taxonomy_*
# 

# Define a default target that collects all the targets from the workflow plus modules
rule all:
    input:
        expand("{project_directory}/data/output/db/{sample_name}.{reference_database}.bestqueryDB", sample_name = config["sample_name"], reference_database = config["reference_database"], project_directory = config["project_directory"]),
        expand("{project_directory}/data/output/db/{sample_name}_queryDB", sample_name = config["sample_name"], project_directory = config["project_directory"]),
        expand("{project_directory}/data/output/{sample_name}.{reference_database}_query_hits.txt", sample_name = config["sample_name"], reference_database = config["reference_database"], project_directory = config["project_directory"]),
        expand("{project_directory}/data/output/{sample_name}.{reference_database}.tax_query_hits.txt", sample_name = config["sample_name"], reference_database = config["reference_database"], project_directory = config["project_directory"]),
        expand("{project_directory}/data/output/{sample_name}.{reference_database}_query_hits.fna", sample_name = config["sample_name"], reference_database = config["reference_database"], project_directory = config["project_directory"])
    wildcard_constraints:
        sample_name = '^[^\/\-]*$',
        reference_database = '^[A-Za-z-]*$'
    default_target: True

wildcard_constraints:
        sample_name = '^[^\/\-]*$',
        reference_database = '^[A-Za-z-]*$'

# This will be for full scale deployment
rule filter_reads:
    input:
        config["sample_directory"] + "/{sample_name}.fastq.gz"
    output:
        "{project_directory}/data/output/{sample_name}.filt.fq.gz"
    threads: 8
    resources:
        mem_mb = 82000
    shell:
        """
        /home/glbrc.org/millican/repos/trait-mapper/workflow/scripts/filter_reads.sh {input} {output}
        """

rule ecc_reads:
    input:
        "{project_directory}/data/output/{sample_name}.filt.fq.gz"
    output:
        "{project_directory}/data/output/{sample_name}.ecc.fq.gz"
    threads: 8
    resources:
        mem_mb = 20000
    shell:
        """
        /home/glbrc.org/millican/repos/trait-mapper/workflow/scripts/ecc_reads.sh {input} {output}
        """

rule query_DB:
    input:
        "{project_directory}/data/output/{sample_name}.ecc.fq.gz"
    output:
        "{project_directory}/data/output/db/{sample_name}_queryDB"
    threads: 1
    resources:
        mem_mb = 24000
    shell:
        """
        /home/glbrc.org/millican/repos/trait-mapper/workflow/scripts/create_db.sh {input} {output}
        """
