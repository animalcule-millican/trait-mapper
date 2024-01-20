#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh
export MMSEQS_FORCE_MERGE=1

mmseqs createdb $1 $2
mmseqs createindex $1 $TMPDIR/tmp --search-type 2 --translation-table 11 --threads $2 --remove-tmp-files 1
