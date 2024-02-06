#!/bin/bash
source /mnt/bigdata/linuxhome/millican/.bashrc
mamba activate sourmash
source random_directory.sh

sourmash tax metagenome --gather $1 --taxonomy $2 --keep-full-identifiers > $3