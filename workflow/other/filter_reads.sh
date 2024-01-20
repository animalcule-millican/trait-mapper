#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

rqcfilter2.sh -Xmx80g in=$1 out=$2 \
chastityfilter=f \
barcodefilter=f \
tmpdir=$TMPDIR \
scafstats=$TMPDIR/scaffoldStats.txt \
kmerstats=$TMPDIR/kmerStats.txt \
log=$TMPDIR/status.log \
filelist=$TMPDIR/file-list.txt \
stats=$TMPDIR/filterStats.txt \
ihist=$TMPDIR/ihist_merge.txt \
outribo=$TMPDIR/ribo.fq.gz \
reproduceName=$TMPDIR/reproduce.sh \
jni=t \
mapk=15 \
maxns=3 \
maq=3 \
minlen=51 \
rqcfilterdata=$RQCFilterData \
qtrim=r \
trimq=5 \
minlength=40 \
removehuman=f \
keephuman=f \
removedog=f \
removecat=f \
removemouse=f \
aggressivehuman=f \
aggressivemicrobe=f \
aggressive=f \
detectmicrobes=f \
removemicrobes=f \
filterpolya=f \
filterpolyg=0 \
phix=t \
lambda=f \
pjet=t \
sip=f \
chloromap=f \
mitomap=f \
ribomap=f \
removeribo=f \
clumpify=f \
dedupe=f \
sketch=f


rm -r $TMPDIR
mamba deactivate


