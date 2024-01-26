#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh
DIR=/home/glbrc.org/millican/repos/trait-mapper/test/time-testing/cluster
INPUT=/home/glbrc.org/millican/repos/trait-mapper/test/data/${1}

CPU=${2}
START=$(date +%s)
# do something
# start your script work here
mmseqs createdb $1 $TMPDIR/seqDB -v 3
mmseqs cluster $TMPDIR/seqDB $TMPDIR/clusterDB tmp --min-seq-id 0.999 --db-load-mode 2 --threads $CPU --remove-tmp-files 1 -v 3
mmseqs createsubdb $TMPDIR/clusterDB $TMPDIR/seqDB $TMPDIR/repDB --id-mode 1 -v 3
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "cluster,${1},$CPU,$DIFF" >> $DIR/cluster_time_testing.csv

rm -r $TMPDIR