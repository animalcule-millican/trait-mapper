#!/bin/bash
source /mnt/bigdata/linuxhome/millican/.bashrc
mamba activate sourmash
source random_directory.sh

sourmash gather $1 $2 --threshold-bp 1000 --output $3