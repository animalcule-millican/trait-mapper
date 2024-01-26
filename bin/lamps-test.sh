#!/bin/bash
source ~/.bashrc

snakemake -s /home/glbrc.org/millican/repos/trait-mapper/workflow/lamps-test.smk --cores all --jobs 1000 --profile glbrc