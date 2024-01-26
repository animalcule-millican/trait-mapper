#!/bin/bash
source ~/.bashrc
snakemake -s /home/glbrc.org/millican/repos/trait-mapper/workflow/Snakefile --cores all --profile HTCondor