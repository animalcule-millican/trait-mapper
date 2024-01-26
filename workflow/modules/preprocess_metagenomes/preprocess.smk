
# moved to processing just the jgi filtered reads. This step is no needed, but I am leaving it here in case we change back to filtering raw reads.
rule filter_reads:
    input:
        config["sample_directory"] + "/{sample_name}.filter-METAGENOME.fastq.gz"
    output:
        "{data_directory}/reads/{sample_name}.filt.fq.gz"
    params:
        RQCFilterData = config["RQCFilterData"],
        tmpdir = create_tmpdir()
    threads: 16
    resources:
        mem_mb = 30000
    conda:
        "trait-mapper"
    shell:
        """
        cd {params.tmpdir}

        rqcfilter2.sh -Xmx30g in={input} out={output} \
        rqcfilterdata={params.RQCFilterData} \
        barcodefilter=f \
        qtrim=r \
        trimq=5 \
        silvalocal=f \
        usetmpdir=t \
        tmpdir={params.tmpdir} \
        merge=f

        cd ~
        rm -r {params.tmpdir}
        """

rule ecc_reads:
    input:
        config["sample_directory"] + "/{sample_name}.filter-METAGENOME.fastq.gz"
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
        "{data_directory}/reads/{sample_name}.ecc.fq.gz"
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

rule sketch_metagenomes:
    input:
        "{data_directory}/reads/{sample_name}.ecc.fq.gz"
    output:
        "{data_directory}/sample_sketches/{sample_name}.sig.gz"
    threads: 1
    resources:
        mem_mb = 30000
    conda:
        "branchwater"
    shell:
        """
        sourmash sketch dna -p k=21k=31,k=51,scaled=1000,abund --name {wildcards.sample_name} -o {output} {input}
        """