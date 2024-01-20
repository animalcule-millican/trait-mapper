#!/usr/bin/env snakemake -s
from snakemake.utils import min_version
configfile: "/home/glbrc.org/millican/repos/trait-mapper/test/config.yml"

workdir: "/home/glbrc.org/millican/repos/trait-mapper/workflow"

project_directory = config["project_directory"]

rule all:
    input:
        expand("{project_directory}/test/data/{sample_name}.filt.fq.gz", project_directory = config["project_directory"], sample_name = config["sample_name"])


rule filter_reads:
    input:
        config["sample_directory"] + "/{sample_name}.fastq.gz"
    output:
        "{wildcard.project_directory}/test/data/{sample_name}.filt.fq.gz"
    params:
        script = config["script_directory"] + "/filter_reads.sh"
    threads: 8
    resources:
        mem_mb = 82000
    script:
        """
        /home/glbrc.org/millican/repos/trait-mapper/workflow/scripts/submit-filter.py {params.script} {input} {output} {threads} {resources.mem_mb}
        """