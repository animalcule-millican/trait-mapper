#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh
export time_table=/home/glbrc.org/millican/projects/trait-mapper/search_time_test_outputs.txt
export targetDB=/home/glbrc.org/millican/ref_db/pgp/nitrogen-cycle/nitrogen-cycle_DB
export sample_name=$1
export CPU=$2
export fasta=/home/glbrc.org/millican/projects/pgp/test_inputs/${sample_name}.fastq.gz
export queryDB=/home/glbrc.org/millican/projects/trait-mapper-test-v3/data/query_db/${sample_name}.queryDB
export repDB=/home/glbrc.org/millican/projects/trait-mapper-test-v3/data/query_db/${sample_name}.repDB
declare -a arr1=(10 100 1000)
# If query database does not exist, create it
if [ ! -f $queryDB ]; then
    echo "No query database found"
    exit 1
fi
mmseqs touchdb $targetDB --threads $CPU -v 3
################################################################################################################################################
for i in "${arr1[@]}"; do
    rm -rf $TMPDIR/*
    max_seq=$i
    START=$(date +%s)
    # do something
    # start your script work here
    mmseqs search $queryDB $targetDB $TMPDIR/resultDB $TMPDIR/tmp --start-sens 1 --sens-steps 3 -s 7 -e 0.1 --max-seqs $max_seq --db-load-mode 2 --translation-table 11 --threads $CPU -v 3
    END=$(date +%s)
    DIFF=$(( $END - $START ))
    echo "search,${sample_name},$CPU,$max_seq,$DIFF" >> $time_table
done
rm -rf $TMPDIR
