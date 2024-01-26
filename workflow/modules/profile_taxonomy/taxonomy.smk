

rule fastgather:
    input:
        sample = "{data_directory}/sample_sketches/{sample_name}.sig.gz",
        ref = "{data_directory}/sketch_list/{taxa}_sketch_list.txt"
    output:
        "{data_directory}/gather/{sample_name}_{taxa}_k{kmer}.fastgather.csv"
    threads: 16
    resources:
        mem_mb = 30000
    conda:
        "branchwater"
    shell:
        """
        sourmash scripts fastgather -o {output} -t 10000 -k {wildcards.kmer} -c {threads} {input.sample} {input.ref}
        """

rule gather:
    input:
        sample = "{data_directory}/sample_sketches/{sample_name}.sig.gz",
        ref = "{data_directory}/sketch/{taxa}_sigs.zip",
        gather = "{data_directory}/gather/{sample_name}_{taxa}_k{kmer}.fastgather.csv"
    output:
        "{data_directory}/gather/{sample_name}_{taxa}_k{kmer}.gather.csv"
    threads: 1
    resources:
        mem_mb = 20000
    conda:
        "branchwater"
    shell:
        """
        sourmash gather {input.sample} {input.ref} --picklist {input.gather}:match_name:ident -o {output} -k {wildcards.kmer} --threshold-bp 10000
        """

rule taxonomy:
    input:
        gather = "{data_directory}/gather/{sample_name}_{taxa}_k{kmer}.gather.csv",
        lin = config["lineage_database"]
    output:
        "{data_directory}/sample_taxonomy/{sample_name}_{taxa}_k{kmer}.csv"
    threads: 1
    resources:
        mem_mb = 20480
    conda:
        "branchwater"
    shell:
        """
        sourmash tax metagenome --gather {input.gather} --taxonomy {input.lin} --keep-full-identifiers > {output}
        """
