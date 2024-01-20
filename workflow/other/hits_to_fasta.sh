#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

/home/glbrc.org/millican/repos/trait-mapper/workflow/scripts/hits_to_fasta.py $1 $2 

rm -r $TMPDIR