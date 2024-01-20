#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

count_mapped_reads.py -i $1 -o $2 

if [ -f $2 ]; then
    rm $1
fi

rm -r $TMPDIR