
rule search_query:
    input:
        query = "{data_directory}/database/{sample_name}.queryDB",
        target = config["reference_directory"] + "/database/pgp/{target_database}_db"
    output:
        hits = "{data_directory}/hit_files/{sample_name}.{target_database}.query_hits.txt"
    params:
        tmpdir = create_tmpdir()
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
