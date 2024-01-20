#!/bin/bash
source ~/.bashrc
mamba activate sourmash
repo=/home/glbrc.org/millican/repos/Slime_Py
seq=/home/glbrc.org/millican/metagenome/lamps/data/reads/raw/mt5016.fastq.gz
#seq=/home/glbrc.org/millican/metagenome/lamps/data/reads/raw/S5201_3.fastq.gz
sig=$repo/test/test_output/mt5016.sig.gz
idx=/home/glbrc.org/millican/ref_db/sourDB/index/viral/viral_21.sbt.zip
outcsv=$repo/test/test_output/fast_gather.csv
ref=/home/glbrc.org/millican/ref_db/sourDB/viral21.txt
find /home/glbrc.org/millican/ref_db/sourDB/signatures/viral/k21 -name "*.sig.gz" > $ref

sourmash sketch dna $seq -p k=21,k=31,k=51,scaled=1000,abund --output $sig
sourmash scripts fastgather $sig $ref -t 1000 -k 21 -o $outcsv
sourmash gather $sig $idx --threshold-bp 1000 --picklist ${outcsv}:match_name:ident --save-matches $repo/test/test_output/gather.matches.zip --ksize 21 --output $repo/test/test_output/gather.csv
