#!/bin/bash
source /home/glbrc.org/millican/.bashrc

snakemake -s /home/glbrc.org/millican/repos/trait-mapper/workflow/map_traits.smk --profile HTCondor
