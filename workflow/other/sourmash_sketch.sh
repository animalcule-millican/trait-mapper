#!/bin/bash
source /mnt/bigdata/linuxhome/millican/.bashrc
mamba activate sourmash
source random_directory.sh

sourmash sketch dna $1 -p k=21,k=31,k=51,scaled=1000,abund --output $2
