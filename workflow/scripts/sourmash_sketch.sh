#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

sourmash sketch dna $1 -p k=21,k=31,k=51,scaled=1000,abund --output $2

trap handle_signal EXIT
trap handle_signal SIGTERM
trap handle_signal ERR
trap handle_signal SIGINT
trap handle_signal SIGQUIT
trap handle_signal SIGKILL