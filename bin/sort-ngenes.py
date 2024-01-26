#!/usr/bin/env python3
from Bio import SeqIO
import pickle
from multiprocessing import Pool

#n_dict = {}
#with open("/home/glbrc.org/millican/ref_db/map_files/ncycle-gene-map.tsv", 'r') as f:
#    for row in f:
#        line = row.strip()
#        header, gene = line.split('\t')
#        n_dict[header] = gene

#with open("/home/glbrc.org/millican/ref_db/map_files/ncycle-gene-map.pkl", 'wb') as f:
#    pickle.dump(n_dict, f)

with open("/home/glbrc.org/millican/ref_db/map_files/ncycle-gene-map.pkl", 'rb') as f:
    n_dict = pickle.load(f)

def sort_seqs(input_rec):
    if input_rec.id in n_dict.keys():
            gene = n_dict[input_rec.id]
            with open(f"/home/glbrc.org/millican/ref_db/pgp/nitrogen-cycle/{gene}.faa", 'a') as handle:
                SeqIO.write(input_rec, handle, 'fasta')

with open("/home/glbrc.org/millican/ref_db/pgp/nitrogen-cycle/nitrogen-cycle.faa", 'r') as f:
    with Pool(processes = 8) as pool:
        pool.map(sort_seqs, SeqIO.parse(f, 'fasta'))