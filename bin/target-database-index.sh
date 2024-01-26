#!/bin/bash
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env
mamba activate trait-mapper
source random_directory.sh

targetDB=/home/glbrc.org/millican/ref_db/pgp/plant-signal-production/plant-signal-production_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/plant-signal-degradation/plant-signal-degradation_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/plant-colonization/plant-colonization_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/root-colonization/root-colonization_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/spore-production/spore-production_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/abiotic-stress-control/abiotic-stress-control_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/plant-stimulation/plant-stimulation_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/antimicrobial-production/antimicrobial-production_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/osmotic-stress/osmotic-stress_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/nutrient-acquisition/nutrient-acquisition_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/nodulation/nodulation_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/microbial-stress-responses/microbial-stress-responses_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/biocontrol/biocontrol_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/bioremediation/bioremediation_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/herbicidal/herbicidal_DB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

targetDB=/home/glbrc.org/millican/ref_db/pgp/phosphorus-cycle/phosphorus-cycle_DB
fasta=/home/glbrc.org/millican/ref_db/pgp/phosphorus-cycle/phosphorus-cycle.faa
mmseqs createdb $fasta $targetDB
mmseqs createindex $targetDB $TMPDIR/tmp --threads 16 --remove-tmp-files 1

rm -r $TMPDIR
