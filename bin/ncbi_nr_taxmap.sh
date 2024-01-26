#!/bin/bash
source ~/.bashrc
random_directory.sh
mmseqs=/home/glbrc.org/millican/mambaforge/envs/trait-mapper/bin/mmseqs
db=/home/glbrc.org/millican/ref_db/trait_db/NT/NTDB
taxmap=/home/glbrc.org/millican/ref_db/trait_db/NTDB_tax_mapping
taxdump=/home/glbrc.org/millican/ref_db/taxdmp
awk 'NR>1 {print $6, $1}' /home/glbrc.org/millican/ref_db/trait_db/gene2accession | awk '$1 != "-"' > $taxmap
$mmseqs createtaxdb $db $TMPDIR/tmp --ncbi-tax-dump $taxdump --tax-mapping-file $taxmap --threads 12
rm -r $TMPDIR