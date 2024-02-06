#!/usr/bin/env snakemake -s
from snakemake.utils import min_version


#rule all:
#    input:
#        expand("{project_directory}/data/output/{sample_name}.{reference_database}.tax.query_hits.txt", sample_name = str(config["sample_name"]), reference_database = str(config["reference_database"]), project_directory = str(config["project_directory"]))

rule sourmash_sketch:
    input:
        "{project_directory}/data/output/{sample_name}.ecc.fq.gz"
    output:
        "{project_directory}/data/output/{sample_name}.sig.gz"
    threads: 1
    resources:
        mem_mb = 20480
    shell:
        """
        scripts/sourmash_sketch.sh {input} {output}
        """

rule sourmash_gather:
    input:
        "{project_directory}/data/output/{sample_name}.sig.gz"
    output:
        "{project_directory}/data/output/{sample_name}-{taxa}.gather.csv"
    params:
        sig_database = "{sig_database}"
        taxa = "{taxa}"
    wildcard_constraints:
        sig_database = "|".join(config["signature_database"])
        taxa = "|".join(config["taxa"])
    threads: 1
    resources:
        mem_mb = 20480
    shell:
        """
        scripts/sourmash_gather.sh {input} {output} {params.sig_database}
        """

rule sourmash_taxonomy:
    input:
        "{project_directory}/data/output/{sample_name}-{taxa}.gather.csv"
    output:
        "{project_directory}/data/output/{sample_name}-{taxa}.taxonomy.csv"
    params:
        tax_database = "{tax_database}"
        taxa = "{taxa}"
    wildcard_constraints:
        tax_database = "|".join(config["taxonomy_database"])
        taxa = "|".join(config["taxa"])
    threads: 1
    resources:
        mem_mb = 20480
    shell:
        """
        scripts/sourmash_gather.sh {input} {output} {params.tax_database}
        """