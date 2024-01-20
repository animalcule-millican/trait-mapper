#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

count_mapped_reads.py -i $1 -o $2 

if [ -f $2 ]; then
    rm $1
fi

trap handle_signal EXIT
trap handle_signal SIGTERM
trap handle_signal ERR
trap handle_signal SIGINT
trap handle_signal SIGQUIT
trap handle_signal SIGKILL
