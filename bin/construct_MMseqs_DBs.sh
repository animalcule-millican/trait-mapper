#!/bin/bash
source /home/glbrc.org/millican/projects/trait-mapper/lamps/etc/trait-mapper.env
mamba activate trait-mapper
# Assign variables

NAME=$(basename -s .faa $1)
DB=/mnt/bigdata/linuxhome/millican/ref_db/pgp/$NAME/${NAME}_DB

mmseqs createdb $1 $DB