#!/usr/bin/env python3
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser = argparse.ArgumentParser(description='Converts a tab delimited file to fasta format')
parser.add_argument('-i','--input', help='Input file')
parser.add_argument('-o', '--output', help='Output file')
parser.add_argument('-g', '--gene', help='gene name')
args = parser.parse_args()

def convert_tab_to_fasta(input_file, output_file, gene_name):
    with open(input_file, 'r') as f, open(output_file, 'a') as o:
        for row in f:
            line = row.strip()
            header, seq = line.split('\t')
            header = f"{header}|{gene_name}"
            rec = SeqRecord(Seq(seq), header, '', '')
            SeqIO.write(rec, o, 'fasta')

if __name__ == '__main__':
    convert_tab_to_fasta(args.input, args.output)
