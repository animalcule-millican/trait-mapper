#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

hits_to_fasta.py $1 $2 

rm -r $TMPDIR

trap handle_signal EXIT
trap handle_signal SIGTERM
trap handle_signal ERR
trap handle_signal SIGINT
trap handle_signal SIGQUIT
trap handle_signal SIGKILL
