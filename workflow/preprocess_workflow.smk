#!/usr/bin/env snakemake -s

rule filter_reads:
    input:
        config["sample_directory"] + "/{sample_name}.fastq.gz"
    output:
        "{project_directory}/data/output/{sample_name}.filt.fq.gz"
    threads: 16
    resources:
        mem_mb = 30000
    shell:
        """
        scripts/filter_reads.sh {input} {output}
        """

rule ecc_reads:
    input:
        "{project_directory}/data/output/{sample_name}.filt.fq.gz"
    output:
        "{project_directory}/data/output/{sample_name}.ecc.fq.gz"
    threads: 16
    resources:
        mem_mb = 30000
    shell:
        """
        scripts/ecc_reads.sh {input} {output}
        """

rule query_DB:
    input:
        "{project_directory}/data/output/{sample_name}.ecc.fq.gz"
    output:
        "{project_directory}/data/output/{sample_name}.queryDB"
    threads: 12
    resources:
        mem_mb = 30000
    shell:
        """
        scripts/create_db.sh {input} {output}
        """
