#!/usr/bin/env snakemake -s

rule search_query:
    input:
        query = "{project_directory}/data/output/{sample_name}.queryDB"
    output:
        hits = "{project_directory}/data/output/{sample_name}.{reference_database}.query_hits.txt",
        best = "{project_directory}/data/output/{sample_name}.{reference_database}.bestqueryDB"
    params:
        target = config["reference_directory"] + "/{reference_database}/{reference_database}_DB"
    threads: 32
    resources:
        mem_mb = 60000
    shell:
        """
        scripts/query_search.sh {input.query} {params.target} {output.best} {output.hits} {threads}
        """

rule search_query_fasta:
    input:
        "{project_directory}/data/output/{sample_name}.{reference_database}.query_hits.txt"
    output:
        "{project_directory}/data/output/{sample_name}.{reference_database}.query_hits.fna"
    threads: 1
    resources:
        mem_mb = 6144
    shell:
        """
        scripts/hits_to_fasta.sh {input} {output}
        """
