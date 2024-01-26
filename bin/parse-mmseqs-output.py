#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

def parse_mmseq_to_fasta(input_file, output_file):
    with open(input_file, 'r') as f, open(output_file, 'a') as o:

        i = 0
        for i,row in enumerate(f):
            i += 1
            if i == 1:
                continue
            line = row.strip()
            line = line.split('\t')
            head_split = line[0].split()
            if head_split[2].startswith('n'):
                seq_name = head_split[1]
            else:
                seq_name = f"{head_split[1]} {head_split[2]}"
            seq_desc = f"{line[2]}|{line[3]}"
            record = SeqRecord(Seq(line[1]), id=head_split[0], name=seq_name, description=seq_desc)
            SeqIO.write(record, o, 'fasta')

def main():
    parser = argparse.ArgumentParser(description='Convert mmseqs2 output to fasta and tabular format')
    parser.add_argument('-i', '--input', help='mmseqs2 output file', required=True)
    parser.add_argument('-o', '--output', help='output fasta file', required=True)
    args = parser.parse_args()
    parse_mmseq_to_fasta(args.input, args.output)