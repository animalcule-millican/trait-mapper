#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

bbcms.sh in=$1 out=$2 bits=4 hashes=3 k=31

if [ -f $2 ]; then
    rm $1
fi