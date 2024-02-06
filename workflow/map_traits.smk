import os
import glob
import random
import pickle

configfile: "/home/glbrc.org/millican/repos/trait-mapper/workflow/config.yml"
workdir: config["working_directory"]

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

def get_sample_names(config):
    refdir = config["sample_directory"]
    g_files = []
    for file in glob.glob(f"{refdir}/*.filter-METAGENOME.fastq.gz"):
        name_file = os.path.basename(file).replace(".filter-METAGENOME.fastq.gz", "")
        g_files.append(name_file)
    return g_files

#g_files = get_reference_names(config)
sample_names = get_sample_names(config)

wildcard_constraints:
    sample_name = "|".join(sample_names),
    data_directory = config["data_directory"],
    sample_directory = config["sample_directory"],
    target_database = "|".join(config["target_database"]),
    reference_directory = config["reference_directory"],
    taxa = "|".join(config["taxa"]),
    kmer = "|".join(config["kmer"]),
    lineage_database = config["lineage_database"]

rule all:
    input:
        expand("{data_directory}/count/{sample_name}.csv", sample_name = sample_names, data_directory = config["data_directory"]),
        #expand("{data_directory}/database/{sample_name}.queryDB", sample_name = sample_names, data_directory = config["data_directory"]),
        #expand("{data_directory}/sample_taxonomy/{sample_name}_{taxa}_k{kmer}.csv", sample_name=sample_names, taxa=config["taxa"], kmer=config["kmer"], data_directory=config["data_directory"]),
        expand("{data_directory}/hit_files/{sample_name}.{target_database}.query_hits.txt", sample_name = sample_names, target_database = config["target_database"], data_directory = config["data_directory"])
    default_target: True


rule ecc_reads:
    input:
        expand("{sample_directory}/{{sample_name}}.filter-METAGENOME.fastq.gz", sample_directory = config["sample_directory"])
    output:
        "{data_directory}/reads/{sample_name}.ecc.fq.gz"
    params:
        tmpdir = create_tmpdir()
    threads: 16
    resources:
        mem_mb = 30000
    conda: 
        "trait-mapper"
    shell:
        """
        export TMPDIR={params.tmpdir}
        bbcms.sh in={input} out={output} bits=4 hashes=3 k=31
        rm -r {params.tmpdir}
        """

rule query_DB:
    input:
        expand("{sample_directory}/{{sample_name}}.filter-METAGENOME.fastq.gz", sample_directory = config["sample_directory"])
    output:
        "{data_directory}/database/{sample_name}.queryDB"
    params:
        tmpdir = create_tmpdir()
    threads: 12
    resources:
        mem_mb = 30000
    conda:
        "trait-mapper"
    shell:
        """
        export TMPDIR={params.tmpdir}
        mmseqs createdb {input} {output}
        mmseqs createindex {output} $TMPDIR/tmp --search-type 2 --translation-table 11 --threads {threads} --remove-tmp-files 1
        rm -r {params.tmpdir}
        """

rule search_query:
    input:
        query = "{data_directory}/database/{sample_name}.queryDB",
        target = config["reference_directory"] + "/{target_database}_db"
    output:
        hits = "{data_directory}/hit_files/{sample_name}.{target_database}.query_hits.txt"
    params:
        tmpdir = create_tmpdir()
    threads: 24
    resources:
        mem_mb = 200000
    conda:
        "trait-mapper"
    shell:
        """
        export TMPDIR={params.tmpdir}
        mmseqs search {input.query} {input.target} $TMPDIR/result $TMPDIR/tmp -e 1.000E-04 -s 3.5 --db-load-mode 3 --alignment-mode 1 --exact-kmer-matching 1 --max-accept 5 --max-rejected 5 --max-seqs 200 --threads {threads} --remove-tmp-files 1 --translation-table 11
        mmseqs filterdb $TMPDIR/result $TMPDIR/besthits --extract-lines 1 --threads {threads}
        mmseqs convertalis {input.query} {input.target} $TMPDIR/besthits {output.hits} --format-mode 0 --db-load-mode 3 --format-output query,qheader,qseq,target,theader,tseq,pident,evalue --threads {threads}
        rm -rf {params.tmpdir}
        """

rule combine_hits:
    input:
        expand("{data_directory}/hit_files/{{sample_name}}.{target_database}.query_hits.txt", target_database = config["target_database"], data_directory = config["data_directory"])
    output:
        "{data_directory}/hit_files/{sample_name}.csv"
    params:
        map_pickle = config['pgp_map_pickle']
    threads: 1
    resources:
        mem_mb = 10000
    conda:
        "trait-mapper"
    shell:
        """
        scripts/combine_hits.py {input} {params.map_pickle} {output}
        """

rule tally_counts:
    input:
        "{data_directory}/hit_files/{sample_name}.csv"
    output:
        "{data_directory}/count/{sample_name}.csv"
    threads: 1
    resources:
        mem_mb = 10000
    conda:
        "trait-mapper"
    shell:
        """
        Rscript scripts/tally_hits.r {input} {output}
        """