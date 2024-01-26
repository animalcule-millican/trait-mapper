#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh
export sample_name=$1
export fasta=/home/glbrc.org/millican/projects/pgp/test_inputs/${sample_name}.fastq.gz
export queryDB=/home/glbrc.org/millican/projects/trait-mapper-test-v3/data/query_db/${sample_name}.queryDB
export repDB=/home/glbrc.org/millican/projects/trait-mapper-test-v3/data/query_db/${sample_name}.repDB
# If query database does not exist, create it
if [ ! -f $queryDB ]; then
    mmseqs createdb $fasta $queryDB
fi
if [ ! -f $repDB ]; then
    mmseqs cluster $queryDB $TMPDIR/clusterDB tmp --min-seq-id 0.8
    mmseqs createsubdb $TMPDIR/clusterDB $queryDB $repDB --id-mode 1
    rm -rf $TMPDIR/clusterD*
fi

