#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

search_output=/home/glbrc.org/millican/repos/trait-mapper/nitrogen-cycle/${1}.txt
fasta_output=/home/glbrc.org/millican/repos/trait-mapper/nitrogen-cycle/${1}.faa
target_DB=/home/glbrc.org/millican/ref_db/trait_db/UniRef90/UniRef90DB

mmseqs createdb /home/glbrc.org/millican/ref_db/pgp/nitrogen-cycle/${1}.faa $TMPDIR/queryDB
mmseqs search $TMPDIR/queryDB $target_DB $TMPDIR/resultDB $TMPDIR/tmp --start-sens 2 -s 7.5 --sens-steps 3
mmseqs filterdb $TMPDIR/resultDB $TMPDIR/besthitDB --comparison-operator le --comparison-value 0.0001 --filter-column 4
mmseqs convertalis $TMPDIR/queryDB $target_DB $TMPDIR/besthitDB $search_output --format-mode 0 --format-output qheader,qseq,theader,tseq

/home/glbrc.org/millican/repos/trait-mapper/bin/tab2fasta.py -i $search_output -o $fasta_output -g $1

rm -r $TMPDIR
#/home/glbrc.org/millican/repos/trait-mapper/bin/expand_ncycle_genes.sh