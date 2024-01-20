#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh
# input definitions
## $1 = input fasta file
## $2 = output reads database for mmseqs
mmseqs createdb $1 $2
mmseqs createindex $2 $TMPDIR/tmp --search-type 2 --translation-table 11 --remove-tmp-files 1
rm -r $TMPDIR