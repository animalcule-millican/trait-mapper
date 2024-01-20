#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

export MMSEQS_FORCE_MERGE=1
query=$1 # query sequence database
target=$2 # besthits database
maptxt=$3 # final output file from mapped reads
cpu=$4 # number of threads

mmseqs search $query $target $TMPDIR/resultsDB $TMPDIR/tmp --search-type 3 -s 7.5 --start-sens 1 --sens-steps 2 --db-load-mode 3 --threads $cpu --remove-tmp-files 1
mmseqs convertalis $query $target $TMPDIR/resultsDB $maptxt --format-mode 0 --db-load-mode 3 --format-output query,qheader,target,theader --threads $cpu

mv $target $TMPDIR
mv ${target}* $TMPDIR
rm -r $TMPDIR