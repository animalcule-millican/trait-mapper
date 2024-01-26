#!/bin/bash
source /home/glbrc.org/millican/.bashrc
mamba activate trait-mapper
source random_directory.sh
export home=/home/glbrc.org/millican
ncyc=/home/glbrc.org/millican/ref_db/pgp/nitrogen-cycle/nitrogen-cycle_DB
out=/home/glbrc.org/millican/ref_db/pgp/expand/nitrogen_cycle
targetDB=/mnt/bigdata/linuxhome/millican/ref_db/trait_db/UniRef100/UniRef100DB

mmseqs search $ncyc $targetDB $TMPDIR/resultDB $TMPDIR/tmp --remove-tmp-files --max-seqs 10000
mmseqs filterdb $TMPDIR/resultDB $TMPDIR/besthitDB --comparison-operator le --comparison-value 0.0001 --filter-column 4
mmseqs convertalis $ncyc $targetDB $TMPDIR/besthitDB $out/besthit.txt --format-mode 4 --format-output theader,tseq,taxid,taxlineage,pident,evalue

/home/glbrc.org/millican/repos/trait-mapper/bin/parse-mmseqs-output.py -i $out/besthit.txt -o $out/expanded_nitrogen_seqs.faa

rm -r $TMPDIR
