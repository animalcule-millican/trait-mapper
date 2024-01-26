#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh
# check input variable if equal to "UniProtKB/Swiss-Prot" then set DB_NAME to "SwissProt"
if [[ $1 == "UniProtKB/Swiss-Prot" ]]; then
    DB_NAME="SwissProt"
# else set DB_NAME to $1
else
    DB_NAME=$1
fi

# set output directory
OUT=/home/glbrc.org/millican/ref_db/trait_db/$DB_NAME
# check if output directory exists and create if not
if [[ ! -d $OUT ]]; then
    mkdir -p $OUT
fi
# run mmseqs databases and download database and create DB
mmseqs databases $1 $OUT/${DB_NAME}DB $TMPDIR/tmp --remove-tmp-files 1 --compressed 1

rm -r $TMPDIR