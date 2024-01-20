#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

export MMSEQS_FORCE_MERGE=1

query=$1
target=$2
repDB=$3
besttxt=$4
cpu=$5

mmseqs search $query $target $TMPDIR/result $TMPDIR/tmp --start-sens 1 --sens-steps 3 -s 7 --db-load-mode 3 --merge-query 1 --threads $cpu --remove-tmp-files 1
mmseqs filterdb $TMPDIR/result $TMPDIR/bestDB --extract-lines 1 --threads $cpu
mmseqs convertalis $query $target $TMPDIR/bestDB $besttxt --format-mode 0 --db-load-mode 3 --format-output query,qheader,qseq,target,theader,tseq,pident,evalue,tcov,qlen,tlen,qset,qsetid,tset,tsetid --threads $cpu
mmseqs createsubdb $TMPDIR/bestDB $query $repDB
mmseqs createindex $repDB $TMPDIR/tmp --search-type 2 --translation-table 11 --remove-tmp-files 1

rm -r $TMPDIR