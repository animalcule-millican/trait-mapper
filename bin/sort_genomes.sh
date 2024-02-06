#!/bin/bash
source ~/.bashrc

export dir=/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_genomes/$1

rename_file()
{
    basename=$(basename "$1")
    # Check if the file name contains the pattern ".genomic.fna.gz"
    if [[ $basename == *".genomic.fna.gz"* ]]; then
        # Replace the pattern ".genomic.fna.gz" with "_genomic.fna.gz"
        newname="${basename//.genomic.fna.gz/_genomic.fna.gz}"
        # Rename the file
        mv "$1" "$dir/$newname"
    fi
}

export -f rename_file

find $dir -maxdepth 1 -type f -name "*.fna.gz" | parallel -j 8 rename_file {}
