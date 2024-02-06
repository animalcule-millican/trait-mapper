#!/usr/bin/env snakemake -s
import os
import glob
import random
import pickle

def create_tmpdir():
    with open("/home/glbrc.org/millican/repos/metagenome_snakemake/etc/adj-aml.pkl", 'rb') as f:
        adj, aml = pickle.load(f)
    temp_dir_base = "/home/glbrc.org/millican/TMPDIR"    # Replace with the base path for temporary directories
    # Construct the temporary directory path
    tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    # Check if the directory exists, and find a new combination if it does
    while os.path.exists(tmpdir):
        tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    # Once we find a combination that does not already exist
    # Create the temporary directory
    os.makedirs(tmpdir, exist_ok=True)
    return tmpdir

wildcard_constraints:
    taxa = "|".join(config["taxa"]),
    database = "|".join(config["database"]),
    refdir = config["reference_directory"],
    index = "|".join(config["index"]),
    genome_index = "|".join(config["genome_file_index"]),
    ref_path = "|".join(config["reference_file_path"]),
    ref_name = "|".join(ref_names),
    target_database = "|".join(config["target_database"])
    #annotation_database = "|".join(config["annotation_database"])

rule search_query:
    input:
        query = "{refdir}/database/genome_file.{genome_index}_db",
        target = "{refdir}/database/pgp/{target_database}_db"
    output:
        hits = "{refdir}/data/output/{sample_name}.{reference_database}.query_hits.txt",
        best = "{refdir}/data/output/{sample_name}.{reference_database}.bestqueryDB"
    threads: 24
    resources:
        mem_mb = 200000
    shell:
        """
        export TMPDIR={params.tmpdir}
        mmseqs search {input.query} {input.target} $TMPDIR/result $TMPDIR/tmp -e 1.000E-04 -s 3.5 --db-load-mode 3 --alignment-mode 1 --exact-kmer-matching 1 --max-accept 5 --max-rejected 5 --max-seqs 200 --threads {threads} --remove-tmp-files 1 --translation-table 11
        mmseqs filterdb $TMPDIR/result {output.best} --extract-lines 1 --threads {threads}
        mmseqs convertalis {input.query} {input.target} {output.best} {output.hits} --format-mode 0 --db-load-mode 3 --format-output query,qheader,qseq,target,theader,tseq,pident,evalue --threads {threads}
        rm -rf {params.tmpdir}
        """

rule search_query_fasta:
    input:
        "{refdir}/data/output/{sample_name}.{reference_database}.query_hits.txt"
    output:
        "{refdir}/data/output/{sample_name}.{reference_database}.query_hits.fna"
    threads: 1
    resources:
        mem_mb = 6144
    shell:
        """
        scripts/hits_to_fasta.sh {input} {output}
        """
