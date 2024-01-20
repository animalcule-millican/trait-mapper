#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

tax=$SOURDB/taxonomy/$2/genbank_taxonomy.csv
sourmash tax metagenome --gather $1 --taxonomy $tax --keep-full-identifiers > $3

trap handle_signal EXIT
trap handle_signal SIGTERM
trap handle_signal ERR
trap handle_signal SIGINT
trap handle_signal SIGQUIT
trap handle_signal SIGKILL