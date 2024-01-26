#!/usr/bin/env python3
from Bio import SeqIO
import pickle

seq_dict = {}
with open("/home/glbrc.org/millican/ref_db/pgp/mgPGPT-db.fasta", 'r') as f:
    for rec in SeqIO.parse(f, 'fasta'):
        seq_dict[rec.id] = {"pgpt": rec.id.split("_")[0], 'seq': rec.seq}

with open("/home/glbrc.org/millican/ref_db/trait-files/pickles/seq_database.pkl", 'wb') as f:
    pickle.dump(seq_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


