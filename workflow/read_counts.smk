#!/usr/bin/env snakemake -s
from snakemake.utils import min_version

rule map_reads:
    input:
        query = "{project_directory}/data/output/{sample_name}.queryDB",
        target = "{project_directory}/data/output/{sample_name}.{reference_database}.bestqueryDB"
    output:
        "{project_directory}/data/output/{sample_name}.{reference_database}.reads_mapped_hits.txt"
    threads: 32
    resources:
        mem_mb = 40000
    shell:
        """
        scripts/map_reads.sh {input.query} {input.target} {output} {threads}
        """
    

rule count_mapped:
    input:
         "{project_directory}/data/output/{sample_name}.{reference_database}.reads_mapped_hits.txt"
    output:
         "{project_directory}/data/output/{sample_name}.{reference_database}.counted_mapped_hits.txt"
    threads: 1
    resources:
        mem_mb = 6000
    shell:
        """
        scripts/count_mapped_reads.sh {input} {output}
        """