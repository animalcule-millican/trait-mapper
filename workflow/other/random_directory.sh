#!/bin/bash
# Define the input file names
source /home/glbrc.org/millican/repos/trait-mapper/etc/trait-mapper.env

# Adjective and animal name files
ADJ=$HOME/repos/trait-mapper/etc/adjectives.txt
AML=$HOME/repos/trait-mapper/etc/names.txt

# Pick two random lines from each file
line_adj="$(shuf -n 1 "$ADJ")"
line_aml="$(shuf -n 1 "$AML")"

# create variable for temporary directory
TMPDIR=$HOME/TMPDIR/${line_adj}-${line_aml}

# checking the random name is not currently in use
if [[ -d $TMPDIR ]]; then # if the directory exists
    while [[ -d $TMPDIR ]]; do # then find new combination until a unique one is found
        line_adj="$(shuf -n 1 "$ADJ")"
        line_aml="$(shuf -n 1 "$AML")"
        # export variable for temporary directory
        export TMPDIR=$HOME/TMPDIR/${line_adj}-${line_aml}
    done
else
    # export variable for temporary directory
    export TMPDIR=$HOME/TMPDIR/${line_adj}-${line_aml}
fi 
# create temporary directory using the randomly assigned names
mkdir -p $TMPDIR